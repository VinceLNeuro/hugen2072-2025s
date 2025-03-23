# Project 4 Work Flow

## 1. Evaluate GWAS hits

### Files:

- GWAS summary stats file: `/ix1/hugen2072-2025s/p4/p4_study1_gwas_results.txt`

- Sample AnnotatedDataFrame: `/ix1/hugen2072-2025s/p4/p4_study1_sample_annotation.RData` -> annot

- Imputed SeqArray GDS file: `/ix1/hugen2072-2025s/p4/p4_study1_imputed_data.gds`

- LocusZoom Plot: `/ix1/hugen2072-2025s/p4/p4_study1_gwas_locuszoom.png`

### Goal of this part: 

Extract the <mark>ref allele dosage</mark> (in this example) for the top hit variant from SeqVarGDS file -> merge it to the phenotype data -> Evalution using sina-plot (`ggforce`)

<br>

### [CheatSheet] Operations on {SeqArray} `SeqVarGDSClass` object

```r
######## 1. Check variant INFO directly ########
variantInfo(gds)

# variant.id chr      pos ref alt
# 1     542766  19 28492903   C   T


######## 2. Set Filters ########

### *Filter by variant id
seqSetFilter(gds, 
             variant.id = my_variant, 
             verbose=TRUE, 
             action = "intersect") # "intersect" – set the current filter to the intersection of selected samples and/or variants

#OUTPUT: # of selected variants: 1


### Filter by Chromosome
seqSetFilterChrom(gds, include=my_chr)


######## 3. Subsetting/Getting Data ########
seqGetData(gds, "field_names")


######## 4. Get dosage ########
SeqVarTools::alleleDosage(gds, n=0) # "n=0 is *** dosage for the reference allele ***"

# same as
SeqVarTools::refDosage(gds)


######## 5. Reset filters ########
seqResetFilter(gds)


######## 6. Close gds file ########
seqClose(gds)
```

<br>

## 2. Aggregate Test

### Workflow
1. Obtain a Null model
    - Same as before!
2. Pick **aggregate type** and create **iterator** for that type
    - sliding window
    - gene or region
3. Run association test
    - Specify **statistical method** ("Burden", "SKAT", "SKATO"…)


### Files:

- Null model is pre-computed in `GENESIS`: `/ix1/hugen2072-2025s/p4/p4_study1_null_model.RData`
    <hr>
- Genotype data (`SeqVarGDSClass`)
- Genotype data + attached Sample Annotation (`SeqVarData` object)


### Analysis for HMGA1 (chr6:34236873-34246231)

- Gene-based Test/Iterator
- All Statistical Methods

```r
#### Gene-based Iterator ####
# Filter to the chromosome of interest
seqSetFilterChrom(seqData, include = "6")
#check
seqGetData(seqData, "chromosome") %>% unique()

# Μake a GRanges object for this gene region +/- 50kb
    # **chr6:34236873-34246231**
gr <- GRanges(seqnames = 6, 
              ranges = IRanges(start=34236873-50e3, end = 34246231+50e3))


#### Statistical Methods ####
# 1. Run a burden test            
iterator <- SeqVarRangeIterator(seqData, gr)
assoc_burden <- assocTestAggregate(iterator, 
                                   null.model, 
                                   AF.max = 1, #not just rare variants
                                   test="Burden")
seqResetFilter(iterator)

assoc_burden$results
assoc_burden$variantInfo[[1]]


# 2. Run a SKAT test
iterator <- SeqVarRangeIterator(seqData, gr)
assoc_SKAT <- assocTestAggregate(iterator, 
                                   null.model, 
                                   AF.max = 1, #not just rare variants
                                   test="SKAT")
seqResetFilter(iterator)

assoc_SKAT$results


# 3. Run a SKAT-O test
iterator <- SeqVarRangeIterator(seqData, gr)
assoc_SKATO <- assocTestAggregate(iterator, 
                                   null.model, 
                                   AF.max = 1, #not just rare variants
                                   test="SKATO")
seqResetFilter(iterator)

assoc_SKATO$results


# Close gds file
seqClose(gds)
```

<br>

## 3. PGS

### Files (pre-harmonized):

- PGS scoring file `/ix1/hugen2072-2025s/p4/p4_pgs002975_scores.txt` 
- PLINK formatted data for study 1 `/ix1/hugen2072-2025s/p4/p4_study1.{bed,bim,fam}`.

```bash
#Q: How many variants matched the PLINK files for the imputed data?

# get variants in pgs file
awk '{print $1}' /ix1/hugen2072-2025s/p4/p4_pgs002975_scores.txt > ls_pgs_variants.txt

# get variants in plink file
awk '{print $2}' /ix1/hugen2072-2025s/p4/p4_study1.bim > ls_plink_bim_variants.txt

# compare the two files 
    # -f: extract pattern from file
    # -F: list of pattern to match
    # -x: match whole line
    # -c: count of matches
grep -Fxf ls_pgs_variants.txt -c ls_plink_bim_variants.txt #182

#### However, the log file says 181 valid predictors
cat apply_pgs-slurm_10238.out | grep -e "^--score.*"
```

### Evalution of the PGS model in this independent data

```r
## Association with height (linear model), adjusting for covariates
eval_lm = lm(height ~ SCORESUM + factor(sex) + factor(pop), data=phenotypes)

summary(eval_lm)

## Report: p-value, partial-R2, adjusted-R2, effect size [95%CI]

```

<br>

## 4. {EasyQC} File Harmonization Before Meta Analysis

`EasyQC` is primarily QC (removing nonsensical data) & checking that alleles are coded consistently between studies

### Files in `/ix1/hugen2072-2025s/p4/`:

- p4_study1_gwas_results.txt
- p4_study2_gwas_results.txt
- **p4_gwas1_gwas2.ecf**
- p4_snp_reference.txt

### EasyQC output major files in `./EasyQC_dir/output/`:

Clean GWAS results
  
  1) `CLEANED.study1.gz`
  2) `CLEANED.study2.gz`

EasyQC report

  3) `p4_gwas1_gwas2.rep`