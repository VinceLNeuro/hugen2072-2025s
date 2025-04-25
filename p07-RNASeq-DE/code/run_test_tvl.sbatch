#!/usr/bin/env bash

######## Slurm resource allocation HUGEN2072 ########
#SBATCH --job-name=rnaseq_alignment
#SBATCH --cluster=teach
#SBATCH --account=hugen2072-2025s
#SBATCH -t 1-00:00 # Runtime in D-HH:MM

#SBATCH --nodes=1 #default - all cores on one machine
#SBATCH --ntasks-per-node=1 #default
#SBATCH --cpus-per-task=2 
#SBATCH --mem=64g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-user=til177@pitt.edu
#SBATCH --mail-type=END,FAIL

#SBATCH --array=0-3 #4 experiments
#SBATCH --output=../output/%x-slurm_%A_%a.out

######## Load software into environment ########
set -ev
module purge
module load cutadapt/2.10
# Load modules for star
module load fastqc/0.11.7
module load star/2.7.5a
module load gcc/9.2.0


######## 0. Setup ArrayJobs ########
exp=$(awk -v lineNumb=$(($SLURM_ARRAY_TASK_ID+1)) 'NR==lineNumb { print $0 }' ls_exp.txt)

## for output file names
exp_prefix=$(echo $exp | sed 's/_ERCC.*//')

## add read{1,2} as input
#       .read1.fastq.gz
#       .read2.fastq.gz
exp_read1="${exp}.read1.fastq.gz"
exp_read2="${exp}.read2.fastq.gz"


echo "Job for Experiment *** ${exp_prefix} *** started at: $(date)"

echo $exp_read1

echo $exp_read2


######## 1. Use cutadapt to trim adapter sequence ########

# Right now the script is set up to take in HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz
# 	and HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
#   and it is set to output in HBR_Rep1.trimmed.R1.fastq.gz and HBR_Rep1.trimmed.R2.fastq.gz
# You need to change the input for each run (directory stays the same)
# You need to change the output file name for each run (use the output file names given in the assignment instructions)
# (You need to change the output file directory from blah/blah/blah/project5/ to a folder of yours; make sure the folder exists)

mkdir -p ../output/1_cutadapt

# -a/-A: Universal Illumina adapters (3')
# -o/-p: OUTPUT (first/second)
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -o ../output/1_cutadapt/${exp_prefix}.trimmed.R1.fastq.gz -p ../output/1_cutadapt/${exp_prefix}.trimmed.R2.fastq.gz \
         ../p7/${exp_read1} ../p7/${exp_read2}


######## 2. Aligning the reads ########

# The inputs are the two files you made in the previous step
# (so remember to update them  each time you run the script)
# (and change the folder from /blah/blah/blah/project5/ of course)
# For the --outFileNamePrefix argument, name the file HBR_Rep1, HBR_Rep2, UHR_Rep1, or UHR_Rep2 depending on which step you're on
# 	also include the path to the folder you're using (i.e., replace /blah/blah/blah/project5)
# Don't change the --genomeDir argument at all

mkdir -p ../output/2_STAR

STAR --runThreadN 2 \
    --runMode alignReads \
    --outSAMstrandField intronMotif \
    --twopassMode Basic \
    --readFilesIn ../output/1_cutadapt/${exp_prefix}.trimmed.R1.fastq.gz \
                  ../output/1_cutadapt/${exp_prefix}.trimmed.R2.fastq.gz \
    --outFileNamePrefix ../output/2_STAR/${exp_prefix} \
    --quantMode GeneCounts \
    --outStd Log \
    --outWigType bedGraph \
    --outWigReferencesPrefix ../output/2_STAR/${exp_prefix}_bedgraph \
    --genomeDir ../p7/STAR_reference/ \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate


#### Main analysis end ####
echo "Job for Experiment *** ${exp_prefix} *** ended at: $(date)  on node $HOSTNAME"

