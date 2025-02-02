# GWASTools GDS Cheat Sheet

## Update Annotation data (`AnnotationDataFrame`, S4)

- `getAnnotation()`/`pData()` (extracts a standard data frame containing the annotation)
- `getMetadata()`/`varMetadata()` (extracts a data frame of variable descriptions)
- `getVariableNames()`/`varLabels()` (extracts a vector of just the variable names)


## Get **Array data** (genotype & intensity{qxy, bl}) (`.gds`) & Combine with Annotation data

- `GdsGenotypeReader("path")` & `GdsIntensityReader("path")` to open the `.gds` files that contain the array data. 
  - class: `GdsGenotypeReader`, `GdsIntensityReader`
  - `@data`, `@snpAnnot`, `@scanAnnot`

- Then, `GenotypeData()` & `IntensityData()` to combine them with the SNP- and sample-annotation 
  - Since the functions we'll be using later require this information to be "packaged" together
  - e.g., `GenotypeData(gds_geno, snpAnnot = snpAnnot, scanAnnot = scanAnnot)`

 
## 