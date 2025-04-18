#!/usr/bin/env bash

######## Slurm resource allocation HUGEN2072 ########
#SBATCH --job-name=p05-WGS-hg38-pipeline-part1
#SBATCH --cluster=htc

#SBATCH --time=1-00:00:00 #24hrs

#SBATCH --nodes=1 #default - all cores on one machine
#SBATCH --ntasks-per-node=1 #default
#SBATCH --cpus-per-task=48  #48 cores

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

######## 1. Read QC ########
mkdir -p output/1_readQC

#### 1.1. checks ####
# check head and tail 
zcat p5/p5_1.fastq.gz | head
zcat p5/p5_2.fastq.gz | head

zcat p5/p5_1.fastq.gz | tail
zcat p5/p5_2.fastq.gz | tail

# check total rows (should divide by 4 for # of reads)
zcat p5/p5_1.fastq.gz | wc -l
zcat p5/p5_2.fastq.gz | wc -l

# check if paired-end (count 1s and 2s)
{
    zcat p5/p5_1.fastq.gz | grep -c "/1"
    zcat p5/p5_2.fastq.gz | grep -c "/2"
}

# check if the reads matches throughout the file
diff -s \
    <(zcat p5/p5_1.fastq.gz | grep "/1" | sed -E 's/\/.*//') \
    <(zcat p5/p5_2.fastq.gz | grep "/2" | sed -E 's/\/.*//')


#### 1.2. fastqc ####
mkdir -p output/1_readQC/1.2_fastqc
fastqc p5/p5_1.fastq.gz -t 48 --outdir=./output/1_readQC/1.2_fastqc
fastqc p5/p5_2.fastq.gz -t 48 --outdir=./output/1_readQC/1.2_fastqc


#### 1.3. cut adapt
mkdir -p output/1_readQC/1.3_PostTrim_Reads

# -j 0 : all cores in this slurm job
# -m 10: drop reads shorter than 10 bp
# -q 20: trim bases off *from the end of the read* until the read has an average quality > 20
# -a/-A: Universal Illumina adapters (3')
# -o/-p: OUTPUT (first/second)
cutadapt -j 0 \
        -m 10 \
        -q 20 p5/p5_1.fastq.gz p5/p5_2.fastq.gz \
        -a AGATCGGAAGAG -A AGATCGGAAGAG \
        -o output/1_readQC/1.3_PostTrim_Reads/1_trimmed.fastq.gz -p output/1_readQC/1.3_PostTrim_Reads/2_trimmed.fastq.gz


#### 1.4 postTrim fastqc ####
mkdir -p output/1_readQC/1.4_PostTrim_QC
fastqc output/1_readQC/1.3_PostTrim_Reads/1_trimmed.fastq.gz -t 48 --outdir=./output/1_readQC/1.4_PostTrim_QC
fastqc output/1_readQC/1.3_PostTrim_Reads/2_trimmed.fastq.gz -t 48 --outdir=./output/1_readQC/1.4_PostTrim_QC



######## See Part 2 for Alignment & Downstream Analyses ########

echo "Job ended at: $(date) on node $HOSTNAME"

