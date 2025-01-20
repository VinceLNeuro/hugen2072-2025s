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


######## Main analysis ########

#### 2. Copies data.vcf to your working directory
cp -v /ix1/hugen2072-2025s/ca/data.vcf ./


#### 3. data.vcf -> data4.bcf.gz (that is sorted and indexed)
# sort (required to be .vcf/.bcf) -> convert (.bcf.gz)     -> index (-c for bcf files)
bcftools sort data.vcf | bcftools view -Ob > data.bcf.gz
bcftools index -c data.bcf.gz

##check
ls data.* #>data.bcf.gz  data.bcf.gz.csi  data.vcf
# bcftools view data.bcf.gz -G | grep -v "##" | head | column -t #start with chr3
bcftools view data.bcf.gz -GH | awk '{print $1}' | sort -s | uniq -c #check all chr -> only 3 & 4


#### 4. Is filtered to include only positions on chromosome 4 -> `data4.bcf.gz`
bcftools view data.bcf.gz -r4 -Ob > data4.bcf.gz

##check
bcftools view data4.bcf.gz -GH | cut -f1 | sort | uniq -c #only chr4


#### (start of PLINK) 5. Then uses PLINK to create a PLINK binary file set version of data4.bcf.gz called data4.{fam,bim,bed};
plink --bcf data4.bcf.gz --make-bed --out data4


#### 6. Then uses PLINK to update sex variable in the data4.fam file, using the sex variable in sex.txt (without copying sex.txt to your own directory),
sex_auxFile="/ix1/hugen2072-2025s/ca/sex.txt"
# "--update-sex expects a file with FIDs and IIDs in the first two columns, 
#   and sex information (1 or M = male, 2 or F = female, 0 = missing) in the (n+2)th column. 
#   If no second parameter is provided, n defaults to 1."

plink --bfile data4 --update-sex $sex_auxFile \
    --make-bed --out data4_updateSex


#### 7,8,9 PLINK: coupled with the phenotypes in phenotype.txt (without copying phenotype.txt to your own directory—and you should use phenotype.txt as an auxiliary file for the next few tasks without altering data4.fam to include phenotype data), to
# - Calculate the allele frequencies of the markers in data4.{fam,bim,bed} in the cases only and
# - Calculate the allele frequencies of the markers in data4.{fam,bim,bed} in the controls only; then
pheno_auxFile="/ix1/hugen2072-2025s/ca/phenotype.txt"
# cut -d" " -f3 $pheno_auxFile | sort | uniq -c
    # 262 1
    # 262 2

## extract case(=2) / control(=1) -> `keep`
    #"--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis."
awk '$3==2 {print $1,$2}' $pheno_auxFile > ./aux_caseKeep.txt
awk '$3==1 {print $1,$2}' $pheno_auxFile > ./aux_controlKeep.txt

## calculate allele-freq for CASE
plink --bfile data4_updateSex --keep aux_caseKeep.txt --freq --out maf_CASE
## for CONTROL
plink --bfile data4_updateSex --keep aux_controlKeep.txt --freq --out maf_CONTROL


#### 10. Performs a GWAS of the phenotype using logistic regression with  * no covariates *.

# "--pheno causes phenotype values to be read from the 3rd column of the specified space- or tab-delimited file, instead of the .fam or .ped file. 
#   The first and second columns of that file must contain family and within-family IDs, respectively."

plink --bfile data4_updateSex --logistic hide-covar beta --ci 0.95 \
    --pheno $pheno_auxFile \
    --out pheno_results
#`pheno_results.assoc.logistic`


######## Main analysis end ########

echo "Job ended at: $(date) on node $HOSTNAME"

crc-job-stats
