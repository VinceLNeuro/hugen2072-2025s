#!/usr/bin/env bash

######## Slurm resource allocation HUGEN2072 ########
#SBATCH --job-name=p05-WGS-hg38-pipeline-part2-24cores
#SBATCH --cluster=teach
#SBATCH --account=hugen2072-2025s
#SBATCH --time=2-00:00:00 #48hrs

#SBATCH --nodes=1 #default - all cores on one machine
#SBATCH --ntasks-per-node=1 #default
#SBATCH --cpus-per-task=24  #24 cores

#SBATCH --mail-user=til177@pitt.edu
#SBATCH --mail-type=END,FAIL

#SBATCH --output=./output/%x-slurm_%A.out

######## Load software into environment ########
set -ev
module purge
module load fastqc/0.11.9
module load cutadapt/2.10
module load gcc/8.2.0 bwa/0.7.17 samtools/1.21
module load gatk/4.5.0.0
module load gcc/8.2.0 bcftools/1.15.1


######## Main analysis ########
echo "Job started at: $(date)"


######## 2. Alignment ########
mkdir -p output/2_alignment

# -R STR	Complete read group header line (user-input --> add to SAM/BAM [important for gatk variant calling]). 
#           ’\t’ can be used in STR and will be converted to a TAB in the output SAM. 
#           The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’.
dir_trimmed_reads="output/1_readQC/1.3_PostTrim_Reads"

bwa mem -t 24 \
    p5/Homo_sapiens_assembly38.fasta \
    ${dir_trimmed_reads}/1_trimmed.fastq.gz ${dir_trimmed_reads}/2_trimmed.fastq.gz \
    -R "@RG\tID:P5\tLB:P5\tSM:P5\tPL:ILLUMINA" \
    | samtools view -bh | samtools sort \
    > output/2_alignment/sorted.bam
# Create index file
samtools index output/2_alignment/sorted.bam

# [submission] first 10 non-header rows of bam 
samtools view output/2_alignment/sorted.bam | head -n10 > luo_submissionB.txt



######## 3. Alignment Quality Control ########
mkdir -p output/3_alignmentQC

#### 3.1. mark duplicates ####
dir_bam="output/2_alignment"
dir_alignmentQC="output/3_alignmentQC"
gatk MarkDuplicatesSpark -I ${dir_bam}/sorted.bam \
                         -O ${dir_alignmentQC}/dupsmarked.bam

#### 3.2. Base QualityScore Recalibration (BQSR) ####
# a) Generate a table of the recalibration statistics
#   -I: BAM/SAM/CRAM file containing reads
#   -R: Reference sequence file
#   --known-sites: One or more _databases_ of known polymorphic sites used to exclude regions around known polymorphisms from analysis.
gatk BaseRecalibrator -I ${dir_alignmentQC}/dupsmarked.bam \
                      -R p5/Homo_sapiens_assembly38.fasta \
                      -O ${dir_alignmentQC}/BQSR.table \
                      --known-sites p5/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                      --known-sites p5/dbsnp_146.hg38.vcf.gz \
                      --java-options "-XX:ParallelGCThreads=24"
# b) Apply the statistics
gatk ApplyBQSR -R p5/Homo_sapiens_assembly38.fasta \
                -I ${dir_alignmentQC}/dupsmarked.bam \
                --bqsr-recal-file ${dir_alignmentQC}/BQSR.table \
                -O ${dir_alignmentQC}/dupsmarked_cleaned.bam \
                --java-options "-XX:ParallelGCThreads=24"

#### 3.3 Alignment Statistics ####
samtools flagstat ${dir_alignmentQC}/dupsmarked_cleaned.bam -@ 24 \
    > ${dir_alignmentQC}/alignment_statistics.out
samtools depth ${dir_alignmentQC}/dupsmarked_cleaned.bam \
    > ${dir_alignmentQC}/depth_statistics.out

# check alignment stats
cat $dir_alignmentQC/alignment_statistics.out

# check depth stats -- 3rd col is counts
head $dir_alignmentQC/depth_statistics.out
#       calculate Avg.Mapped.ReadDepth
sum_perBaseMappedDepth=$(awk '{sum += $3} END{print sum}' $dir_alignmentQC/depth_statistics.out)
n_bases=$(wc -l $dir_alignmentQC/depth_statistics.out)
avg_mappedDepth=$(awk -v sum="$sum_perBaseMappedDepth" -v n="$n_bases" 'BEGIN{ print sum / n }')
echo "-------- Avg.Mapped.ReadDepth = $avg_mappedDepth --------"



######## 4. Genotyping/Variant-Calling ########
dir_genotyping="output/4_genotyping"
mkdir -p $dir_genotyping

gatk HaplotypeCaller -R p5/Homo_sapiens_assembly38.fasta \
                    -I ${dir_alignmentQC}/dupsmarked_cleaned.bam \
                    -O ${dir_genotyping}/genotypes.g.vcf.gz \
                    -ERC GVCF \
                    -OVI \
                    --native-pair-hmm-threads 24

# gvcf -> vcf
gatk GenotypeGVCFs -R p5/Homo_sapiens_assembly38.fasta \
                -V ${dir_genotyping}/genotypes.g.vcf.gz \
                -O ${dir_genotyping}/genotypes.vcf.gz \
                --java-options "-XX:ParallelGCThreads=24"
#               -V genotypes_2.g.vcf.gz \ # if there is a 2nd participant
#               -V genotypes_3.g.vcf.gz \ # if there is a 3rd participant, etc.
# Check vcf file
bcftools view ${dir_genotyping}/genotypes.vcf.gz | grep -v "##" | head | column -t



######## 5. Genotype QC ########
mkdir -p output/5_genotypeQC
dir_genotypeQC="output/5_genotypeQC"

#### 5.1. Create site-information-only (no geno) VCF ####
gatk MakeSitesOnlyVcf -I ${dir_genotyping}/genotypes.vcf.gz \
                      -O ${dir_genotypeQC}/sites_only.vcf.gz

#### 5.2. Recalibration for SNPs and Indels ####
#   -an: The names of the annotations which should used for calculations
# a) INDEL
gatk VariantRecalibrator -mode INDEL \
                        -R p5/Homo_sapiens_assembly38.fasta \
                        -V ${dir_genotypeQC}/sites_only.vcf.gz \
                        -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
                        -resource:mills,known=false,training=true,truth=true,prior=12 \
                                    p5/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                        -resource:dbsnp,known=true,training=false,truth=false,prior=2 \
                                    p5/dbsnp_146.hg38.vcf.gz \
                        -resource:axiomPoly,known=false,training=true,truth=false,prior=10 \
                                    p5/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
                        --java-options "-XX:ParallelGCThreads=24" \
                        -O ${dir_genotypeQC}/indels.recal \
                        --tranches-file ${dir_genotypeQC}/indels.tranches
# b) SNP
gatk VariantRecalibrator -mode SNP \
                        -R p5/Homo_sapiens_assembly38.fasta \
                        -V ${dir_genotypeQC}/sites_only.vcf.gz \
                        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
                        -resource:hapmap,known=false,training=true,truth=true,prior=15 \
                                p5/hapmap_3.3.hg38.vcf.gz \
                        -resource:omni,known=false,training=true,truth=true,prior=12 \
                                p5/1000G_omni2.5.hg38.vcf.gz \
                        -resource:1000G,known=false,training=true,truth=false,prior=10 \
                                p5/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
                        -resource:dbsnp,known=true,training=false,truth=false,prior=7 \
                                p5/dbsnp_146.hg38.vcf.gz \
                        --java-options "-XX:ParallelGCThreads=24" \
                        -O ${dir_genotypeQC}/snps.recal \
                        --tranches-file ${dir_genotypeQC}/snps.tranches
# c) Apply INDEL
gatk ApplyVQSR -mode INDEL \
            -R p5/Homo_sapiens_assembly38.fasta \
            -V ${dir_genotyping}/genotypes.vcf.gz \
            --recal-file ${dir_genotypeQC}/indels.recal \
            --tranches-file ${dir_genotypeQC}/indels.tranches \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            --java-options "-XX:ParallelGCThreads=24" \
            -O ${dir_genotypeQC}/genotypes_indelqc.vcf.gz
# d) Apply SNP
gatk ApplyVQSR -mode SNP \
            -R p5/Homo_sapiens_assembly38.fasta \
            -V ${dir_genotypeQC}/genotypes_indelqc.vcf.gz \
            --recal-file ${dir_genotypeQC}/snps.recal \
            --tranches-file ${dir_genotypeQC}/snps.tranches \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            --java-options "-XX:ParallelGCThreads=24" \
            -O ${dir_genotypeQC}/final.vcf.gz

# [submission] (d) A text file of the first ten   non-header rows of   the output vcf file from step12(d).
bcftools view ${dir_genotypeQC}/final.vcf.gz | grep -v "##" | head -n10 \
    > luo_submissionD_finalVCF.txt



#### 5.3. Get __Genotyping/Variant-Calling Quality Metrics__ ####
gatk CollectVariantCallingMetrics -I ${dir_genotypeQC}/final.vcf.gz \
                                --DBSNP p5/dbsnp_146.hg38.vcf.gz \
                                -O ${dir_genotypeQC}/genotype_metrics
# genotype_metrics.variant_calling_detail_metrics
# genotype_metrics.variant_calling_summary_metrics

#### VCF-wise metrics
cat ${dir_genotypeQC}/genotype_metrics.variant_calling_summary_metrics
#       details see: https://broadinstitute.github.io/picard/picard-metric-definitions.html
awk 'NR>6 { print $1,$2,$3,$6,$7, $8,$9,$12,$13,$14 }' ${dir_genotypeQC}/genotype_metrics.variant_calling_summary_metrics \
    | column -t



######## Main analysis end ########

echo "Job ended at: $(date) on node $HOSTNAME"

