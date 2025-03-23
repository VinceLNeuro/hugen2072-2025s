#!/usr/bin/env bash

######## Slurm resource allocation HUGEN2072 ########
#SBATCH --job-name=apply_pgs
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
module load plink/1.90b6.7
set -ev

echo "Job started at: $(date)"

#### Main analysis ####

# The 'sum' modifier causes sums (Linear Comb. Sum) to be reported instead.
plink --bfile /ix1/hugen2072-2025s/p4/p4_study1 \
    --score /ix1/hugen2072-2025s/p4/p4_pgs002975_scores.txt header sum \
    --out ./study1_scored

#### Main analysis end ####

echo "Job ended at: $(date) on node $HOSTNAME"

