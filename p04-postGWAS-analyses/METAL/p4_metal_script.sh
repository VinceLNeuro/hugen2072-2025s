# Begin script 

# Set the method to use
SCHEME STDERR       # inverse-variance weighted
STDERR SE           # harmonized colname from EasyQC
GENOMICCONTROL ON   # apply genomic control

# Describe and process the results from study 1
MARKER		cptid
ALLELE		EFFECT_ALLELE	OTHER_ALLELE
PVALUE		PVAL
WEIGHT		N
EFFECT		BETA
SEPARATOR  TAB
COLUMNCOUNTING LENIENT
PROCESS     ../EasyQC_dir/output/CLEANED.study1.gz


# Describe and process from study 2
MARKER		cptid
ALLELE		EFFECT_ALLELE	OTHER_ALLELE
PVALUE		PVAL
WEIGHT		N
EFFECT		BETA
SEPARATOR  TAB
COLUMNCOUNTING LENIENT
PROCESS     ../EasyQC_dir/output/CLEANED.study2.gz


# Run the meta-analysis, calculate heterogeneity, and then quit METAL
ANALYZE HETEROGENEITY
QUIT

