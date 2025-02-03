# GWASTools GDS Cheat Sheet

## 1. Update Annotation data (`-AnnotationDataFrame`, S4)

> getter-function (both sides) / setter-function

- `getAnnotation()`/`pData()` (extracts a standard data frame containing the annotation)
- `getMetadata()`/`varMetadata()` (extracts a data frame of variable descriptions)
- `getVariableNames()`/`varLabels()` (extracts a vector of just the variable names)


## 2. Get **Array data** (genotype & intensity{qxy, bl}) (`.gds`) & Combine with Annotation data

- `GdsGenotypeReader("path")` & `GdsIntensityReader("path")` to open the `.gds` files that contain the array data. 
  - *class*: `GdsGenotypeReader`, `GdsIntensityReader`
  - *functions*: `getX()`, `getY()` for intensity variables data

- Then, `GenotypeData()` & `IntensityData()` to **combine** them with the SNP- and sample-annotation to create objects with the data-type required by functions in the following pipeline

  - example: 
    ```r
    GenotypeData(gds_geno, snpAnnot = snpAnnot, scanAnnot = scanAnnot)
    ```

<br>

## {GWASTools} QC Step 1: Missingness

```r
missingGenotypeBySnpSex(genoData, #GenotypeData object
                        scan.exclude = NULL, #vector, `scanID` to be excluded
                        verbose = TRUE)
```

- Output:
  1. `missing.counts`: A matrix (row = per SNP, column = per sex) containing the number of missing SNP counts for males and females, respectively.
    
  2. `scans.per.sex`: A vector, containing the number of males and females, respectively.
    
  3. **[target] `missing.fraction`**:	A vector (names(.) = SNP), containing the **missingness rate for each SNP**, with females excluded for the Y chromosome.

<br>

```r
missingGenotypeByScanChrom(GenotypeData, 
                           snp.exclude = NULL, #vector, `snpID` to be excluded from missing count
                           verbose = TRUE)
```

- Output: 
  1. `missing.counts`: A matrix (rows = scans, cols = unique chromosomes) containing the **number of missing SNP's** for each scan and chromosome.

  2. `snps.per.chr`: A vector, containing the number of non-excluded SNPs for each chromosome.

  3. `missing.fraction`: A vector, containing the **missingness for each scan** over all chromosomes, excluding the Y chromosome for females.

<br>

## {GWASTools} QC Step 2: Reported Sex vs Inferred Sex (possibly sex-chromosome aneuploidy)

**Data Type** to work on: IntensityData (`qxyData`)
    
```r
meanIntensityByScanChrom(intenData, #IntensityData object
                        vars = c("X", "Y"), #char-vector for intensity variables names
                                            #(two allele-specific probes to determine the SNP genotype)
                        snp.exclude = NULL, 
                        verbose = TRUE)
```
- Output (partial):  
    1. **[target]`$mean.intensity`**: matrix (row=scan, col=chromosome), containing mean of summed intensity values (integrated/summed X probe and Y probe values) for each scan for each chromosome
    
    2. `$mean.X`: matrix containing means of intensity values for **just X probe**

<br>

```r
hetByScanChrom(genoData, snp.exclude = NULL, verbose = TRUE)
``` 
- Output: matrix (row=scan, col=chr), containing <u>chromosomes heterozygosity rates</u>, including a column **"A" for all autosomes**.

<br>

## {GWASTools} QC Step 3: B allele frequency (BAF) and log R ratio (LRR) plots

`chromIntensityPlot()`
* [NOTE] Check the manual for more information

* e.g., *chromosome 22 was flagged for scanID 286.*
    ```r
    chromIntensityPlot(intenData = blData, #IntensityData
                       scan.ids = 286, chrom.ids = 22)
    ```