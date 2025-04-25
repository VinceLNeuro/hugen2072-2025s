#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J run_test
#SBATCH --output=run_test.out
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.
#SBATCH --mem=64g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero

# Use cutadapt to trim adapter sequence
# Right now the script is set up to take in HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz
# 	and HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
#   and it is set to output in HBR_Rep1.trimmed.R1.fastq.gz and HBR_Rep1.trimmed.R2.fastq.gz
# You need to change the input for each run (directory stays the same)
# You need to change the output file name for each run (use the output file names given in the assignment instructions)
# (You need to change the output file directory from blah/blah/blah/project5/ to a folder of yours; make sure the folder exists)
module load cutadapt/2.10
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o /blah/blah/blah/project7/HBR_Rep1.trimmed.R1.fastq.gz \
-p /blah/blah/blah/project7/HBR_Rep1.trimmed.R2.fastq.gz \
/ix1/hugen2072-2025s/p7/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz \
/ix1/hugen2072-2025s/p7/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

# Load modules
module load fastqc/0.11.7
module load star/2.7.5a
module load gcc/9.2.0

# Aligning the reads
# The inputs are the two files you made in the previous step
# (so remember to update them  each time you run the script)
# (and change the folder from /blah/blah/blah/project5/ of course)
# For the --outFileNamePrefix argument, name the file HBR_Rep1, HBR_Rep2, UHR_Rep1, or UHR_Rep2 depending on which step you're on
# 	also include the path to the folder you're using (i.e., replace /blah/blah/blah/project5)
# Don't change the --genomeDir argument at all
STAR \
--runMode alignReads \
--outSAMstrandField intronMotif \
--twopassMode Basic \
--readFilesIn /blah/blah/blah/project7/HBR_Rep1.trimmed.R1.fastq.gz \
              /blah/blah/blah/project7/HBR_Rep1.trimmed.R2.fastq.gz \
--outFileNamePrefix /blah/blah/blah/project7/HBR_Rep1 \
--quantMode GeneCounts \
--outStd Log \
--outWigType bedGraph \
--outWigReferencesPrefix ./bedgraph \
--genomeDir /ix1/hugen2072-2025s/p7/STAR_reference/ \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate

