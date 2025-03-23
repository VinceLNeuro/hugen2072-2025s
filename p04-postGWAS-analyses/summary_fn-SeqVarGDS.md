# Project 4 Work Flow

## 1. Evaluate GWAS hits

### Files:

- GWAS summary stats file: `/ix1/hugen2072-2025s/p4/p4_study1_gwas_results.txt`

- Sample AnnotatedDataFrame: `/ix1/hugen2072-2025s/p4/p4_study1_sample_annotation.RData` -> annot

- Imputed SeqArray GDS file: `/ix1/hugen2072-2025s/p4/p4_study1_imputed_data.gds`

- LocusZoom Plot: `/ix1/hugen2072-2025s/p4/p4_study1_gwas_locuszoom.png`

### *Goal of this part: 

Extract the <mark>ref allele dosage</mark> (in this example) for the top hit variant from SeqVarGDS file -> merge it to the phenotype data -> Evalution using sina-plot (`ggforce`)

<br>


### ***Operations on `SeqVarGDSClass` object***

```r
#### Check variant INFO directly ####
variantInfo(gds)

# variant.id chr      pos ref alt
# 1     542766  19 28492903   C   T


#### Set Filters ####

### *Filter by variant id
seqSetFilter(gds, 
             variant.id = my_variant, 
             verbose=TRUE, 
             action = "intersect") # "intersect" – set the current filter to the intersection of selected samples and/or variants

#OUTPUT: # of selected variants: 1


### Filter by Chromosome
seqSetFilterChrom(gds, include=my_chr)


#### 2. Get dosage ####
SeqVarTools::alleleDosage(gds, n=0) # "n=0 is *** dosage for the reference allele ***"

# same as
SeqVarTools::refDosage(gds)


#### 3. reset filters ####
seqResetFilter(gds)

```

