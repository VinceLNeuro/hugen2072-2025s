#!/usr/bin/env bash

######## Slurm resource allocation HUGEN2072 ########
#SBATCH --job-name=coding_assessment
#SBATCH --cluster=teach
#SBATCH --account=hugen2072-2025s
#SBATCH --time=1:00:00 #1hr

#SBATCH --nodes=1 #default - all cores on one machine
#SBATCH --ntasks-per-node=1 #default
#SBATCH --cpus-per-task=2 
#SBATCH --mem=2000MB #2GB
#SBATCH --mail-user=til177@pitt.edu
#SBATCH --mail-type=END,FAIL

#SBATCH --output=./%x-slurm_%A.out

######## Load software into environment ########
module purge
module load plink/1.90_20230116     gcc/8.2.0 bcftools/1.15.1

set -e
# set -euo pipefail # detect errors
set -v            # print all commands

echo "Job started at: $(date)"


#### Main analysis ####

#### 2. Copies data.vcf to your working directory
cp -v /ix1/hugen2072-2025s/ca/data.vcf ./


#### 3. data.vcf -> data4.bcf.gz (that is sorted and indexed)
# sort (required to be .vcf/.bcf) -> convert (.bcf.gz)     -> index (-c for bcf files)
bcftools sort data.vcf | bcftools view -Ob > data4.bcf.gz
bcftools index -c data4.bcf.gz

# check
ls data4.* #>data4.bcf.gz  data4.bcf.gz.csi
bcftools view data4.bcf.gz -G | head -n14 #start with chr3
bcftools view data4.bcf.gz -GH | awk '{print $1}' | sort -s | uniq -c #check all chr (only 3 & 4)


#### 4. Is filtered to include only positions on chromosome 4




#### Main analysis end ####

echo "Job ended at: $(date) on node $HOSTNAME"

crc-job-stats
