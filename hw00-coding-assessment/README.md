# 0. Background

- 4 files for the project: `/ix1/hugen2072-2025s/ca/`

- `crc-interactive --teach --account hugen2072-2025s -t 3:00:00`

- What you must submit to Canvas:

    - (*DONE*) Your slurm script file for steps 1–10: `assessment_slurm.sh`
        - The slurm log file for the script when you ran the script: `coding_assessment-slurm_6367.out`

    - (*DONE*) The RMarkdown file for steps 11–20: `assessment.Rmd`
        - The knitted HTML file produced by the RMarkdown file: `assessment.html`

<br>

# 1. Tasks

1. (*DONE*) Create and run a script on the CRC cluster that

2. (*DONE*) Copies data.vcf to your working directory,

3. (*DONE*) Then uses bcftools to create a data4.bcf.gz file that is sorted and indexed, and

4. (*DONE*) Is filtered to include only positions on **chromosome 4**;

5. (*DONE*) Then uses PLINK to create a PLINK binary file set version of data4.bcf.gz called data4.{fam,bim,bed};

6. (*DONE*) Then uses PLINK to update sex variable in the data4.fam file using the sex variable in sex.txt (without copying sex.txt to your own directory)
    - `/ix1/hugen2072-2025s/ca/`

7. (*DONE*) Then uses PLINK, coupled with the phenotypes in phenotype.txt (without copying phenotype.txt to your own directory—and you should use phenotype.txt as an auxiliary file for the next few tasks without altering data4.fam to include phenotype data), to

8. (*DONE*) Calculate the allele frequencies of the markers in data4.{fam,bim,bed} in the **cases only** and

9. (*DONE*) Calculate the allele frequencies of the markers in data4.{fam,bim,bed} in the **controls only**; then

10. (*DONE*) Performs a GWAS of the phenotype using logistic regression with no covariates.

11. (*DONE*) Then, create an RMarkdown file that knits to an HTML file which

12. (*DONE*) Reads the two allele frequency files and the logistic regression output into R, and

13. (*DONE*) Plots a histogram of the case allele frequencies,

14. (*DONE*) Plots a histogram of the control allele frequencies,

15. (*DONE*) Plots a scatterplot of the case individuals’ allele frequencies (y axis)

16. (*DONE*) Versus the control individuals’ allele frequencies (x axis),

17. (*DONE*) Plots a scatterplot of −log₁₀(p values) (y axis) versus chr4 basepair position (x axis),

18. (*DONE*) Reports the marker and p value of the marker with the lowest p value; then

19. (*DONE*) Loads the snp_annot.RData workspace and merges the p values into the SNP annotation data frame object snp_annot_df in such a way that

20. (*DONE*) The order of the SNPs in snp_annot_df is unchanged.
