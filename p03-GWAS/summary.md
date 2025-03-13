# Height GWAS

## 0. Data Description

### Sample data

* N = 2548
* median heights: 
  * males > females
  * AFR > EAS > SAS > AMR > EUR
    ```
    pop median_heights
    AFR	173.7158			
    EAS	168.8816			
    SAS	168.6429			
    AMR	165.2727			
    EUR	161.2815
    ```
Based on the EDA, it is clear that `sex` and `pop` are associated with heights. They might also affect genetic info (especially population structure), so we should adjust for them in the GWAS. 

### Genotype Data (Post-imputation)

* From TOPMed -> downsampled

* Variant info (`snp.info`) 
  * 647,629 variants in total
    * = 10,009 genotyped + 637,620 imputed variants 
    * 23,747 variants with ALT allele freq = 0 | 1

  * Determine the imputation R2 threshold (<mark>r2 > 0.3</mark> -> n=241,225)
    * Reason: We can see that the amount of low-quality imputed variants stop to drop at around **r2 = 0.3**, so this could be a reasonable threshold while retaining more reasonably high-quality variants.

  * Determine MAF threshold
    * MAF > 5% (alt <mark>af > 0.05 & af < 0.95</mark>) 

* Kept SNP list (<mark>n = 58,513 SNPS to be included in GWAS</mark>)
    ```r
    snps.keep <- snp.info %>% 
        filter( (r2 > 0.3) 
              & (af > 0.05 & af < 0.95) ) %>% pull(snp.id)

    length(snps.keep) #58513 SNPS
    ```

<br>

## 1. PCA

* From the PC scree plot, the final largest drop is at PC4 -> PC5, which informs the <mark>selection of the first 5 PCs</mark>.

* However, from pairwise PC plot (first 8 PCs), we can see that first 4 PCs are clearly explaining cross-ancestry variances, while PC5, 6, 8 are explaining variances within the AFR population and PC7 is explaining that within the EAS population. It is noteworthy that PC6 vs PC8 seems separating AFR population into 7 subpopulations. Thus, it seems reasonable to include the first 8 PCs.
    * stored in `p03-GWAS/PCs_20250313.pdf`

* Further from the Parellel Coordinate Plot, we can see that PC9 is further separating EUR population from others, and PC10 seems another PC explaining AFR sub population structure. There is no interesting things happen after PC10. 

<mark>Thus, we decided to include the first 10 PCs</mark>



