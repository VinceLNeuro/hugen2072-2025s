#!/usr/bin/env bash

######## Slurm resource allocation HUGEN2072 ########
#SBATCH --job-name=CHANGE_THIS
#SBATCH --cluster=teach
#SBATCH --account=hugen2072-2025s
#SBATCH --time=1:00:00 #1hr

#SBATCH --nodes=1 #default - all cores on one machine
#SBATCH --ntasks-per-node=1 #default
#SBATCH --cpus-per-task=2 
#SBATCH --mem=2000MB #2GB
#SBATCH --mail-user=til177@pitt.edu
#SBATCH --mail-type=END,FAIL

#SBATCH --output=../RELATIVE/PATH/%x-slurm_%A.out

######## Load software into environment ########
# set -euo pipefail # detect errors
set -ev
# module purge
# module load plink/1.90_20230116 #gcc/8.2.0 bcftools/1.15.1


#### Main analysis ####
echo "Job started at: $(date)"


#### Main analysis end ####

echo "Job ended at: $(date) on node $HOSTNAME"

