#!/usr/bin/env bash

######## Slurm resource allocation HUGEN2072 ########
#SBATCH --job-name=sv-calling
#SBATCH --cluster=teach
#SBATCH --account=hugen2072-2025s
#SBATCH --time=1:00:00 #1hr

#SBATCH --nodes=1 #default - all cores on one machine
#SBATCH --ntasks-per-node=1 #default
#SBATCH --cpus-per-task=16

#SBATCH --mail-user=til177@pitt.edu
#SBATCH --mail-type=END,FAIL

#SBATCH --output=./log/%x-slurm_%A.out

######## Load software into environment ########
set -ev
module purge
module load gcc/8.2.0 samtools/1.12 bcftools/1.15.1


echo "Job started at: $(date)"

######## 1. Extract chr22 alignments ########
mkdir -p output

#   -T: required whenever writing __CRAM__
samtools view -T p6/GRCh38_full_analysis_set_plus_decoy_hla.fa \
              -bS \
              -h p6/NA12778.final.cram chr22 \
              -o output/NA12778.chr22.bam \
#   index (required by manta)
samtools index output/NA12778.chr22.bam


######## 2. manta ########

#### 2.1. Configure manta ####
python /ihome/crc/install/manta/manta-1.6.0.centos6_x86_64/bin/configManta.py \
    --referenceFasta p6/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --bam output/NA12778.chr22.bam \
    --runDir output/manta_test
#   Make a directory with a python script for downstream workflow

#### 2.2 Run manta workflow ####
python output/manta_test/runWorkflow.py


######## 3. Extract Variants > 1 kb, Convert vcf to bed Format ########
manta_wkdir="./output/manta_test"

# View the header of `diploidSV.vcf.gz`
bcftools view ${manta_wkdir}/results/variants/diploidSV.vcf.gz | grep -v "##" | head

# Conversion Workflow
#       a) Extract the deletions and duplications from vcf
#       b) Clean VCF format: only get first part of INFO (END pos)--use regexp for sep
#       c) Keep ‘events’ (structural variants) length > 1 kb.
#       d) Create a bed file
zcat ${manta_wkdir}/results/variants/diploidSV.vcf.gz \
    | awk '$5 ~ /DEL/ || $5 ~ /DUP/ {print $0}' \
    | awk -F"[\t;]" '{print $1, $2, $8, $5}' OFS="\t" \
    | sed 's/END=//g' \
    | awk '$3-$2 >= 1000 {print}' \
    > output/gt1kb.cnv.bed


######## 4. Convert chr22 BAM to indexed CRAM --> easier for IGV ########
# output as CRAM `-C` require `-T ref.fa`
samtools view -T p6/GRCh38_full_analysis_set_plus_decoy_hla.fa -C \
              -h output/NA12778.chr22.bam \
              -o output/NA12778.chr22.cram

samtools index output/NA12778.chr22.cram
# samtools quickcheck output/NA12778.chr22.cram
# echo $?



######## Main analysis end ########
echo "Job ended at: $(date) on node $HOSTNAME"

