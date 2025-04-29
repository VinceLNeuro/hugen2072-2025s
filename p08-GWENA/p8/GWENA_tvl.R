# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# install.packages("biomaRt")
# BiocManager::install("GWENA")

# setwd("p08-GWENA/p8/") #if run from RStudio
rm(list=ls())
library(GWENA)
library(magrittr) # Not mandatory, we use the pipe `%>%` to ease readability.
threads_to_use <- 2

################################################################################

#### Q1. explain how the next line deals with the format difference from the example file. ####
# Hint 1: what is the purpose of function t that is added to the example file. 
# Hint 2: you can try the next line without function t to see the difference.
UHRBHR_expr = t(read.table("./UHRBHRExpressions.csv", sep=',', header=TRUE, row.names=1))
# ncol(UHRBHR_expr)
# nrow(UHRBHR_expr) 
dim(UHRBHR_expr)
UHRBHR_expr[1:4,1:4]
is_data_expr(UHRBHR_expr)

# #### if not doing transposition ####
# test_expr_mat = read.table("./UHRBHRExpressions.csv", sep=',', header=TRUE, row.names=1)
# head(test_expr_mat, 4)
# dim(test_expr_mat)
# is_data_expr(test_expr_mat)


# read Metadata
UHRBHR_traits = read.table("UHRBHR_traits.txt", header=TRUE, row.names=1)
UHRBHR_traits
unique(UHRBHR_traits$Condition) #2 cond.


#### Q2. explain how the following module selects expressed genes to further anlayze. ####
dim(UHRBHR_expr)[2]
# Selects only functional genes for downstream analysis 
UHRBHR_expr_filtered <- UHRBHR_expr[,colSums(UHRBHR_expr)>1]
dim(UHRBHR_expr_filtered)[2]


#### Q3. explain what the next line does for the given gene names and why this is needed. ####
# Hint: run the subsequent lines without line 29 and 30 to see the difference, especially in calling bio_enrich
genes<-colnames(UHRBHR_expr_filtered)
genes.adj <- substr(genes, 1, 15) #remove the version number
colnames(UHRBHR_expr_filtered) <- genes.adj
# Remaining number of genes
ncol(UHRBHR_expr_filtered)

# Further downsize to 200 genes
UHRBHR_expr_filtered<-UHRBHR_expr_filtered[,1:200]

# Compute *correlation*
net <- build_net(UHRBHR_expr_filtered, cor_func = "spearman", n_threads = threads_to_use)
    ## Spearman is less sensitive to outliers which are frequent in transcriptomics datasets and does not assume normal distribution.

# power selected: 1
net$metadata$power
fit_power_table <- net$metadata$fit_power_table
# find power law fitted R2
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]


# Detect modules (hierarchical clustering)
modules <- detect_modules(UHRBHR_expr_filtered, 
                          net$network, 
                          detailled_result = TRUE,
                          merge_threshold = 0.25)

# Number of modules before merging :
length(unique(modules$modules_premerge))
# Number of modules after merging: 
length(unique(modules$modules))

### Plot 1
# png(filename = "../submission/layout_modules_PrePostMerge.png", width = 4, height = 4, units = "in", res = 500)
layout_mod_merge <- plot_modules_merge(modules_premerge = modules$modules_premerge, 
                                       modules_merged   = modules$modules)
# dev.off()

### Plot 2
# png(filename = "../submission/numbGenes_PerMod.png", width = 6, height = 6, units = "in", res = 500)
ggplot2::ggplot(data.frame(modules$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")
# dev.off()

### Plot 3
# png(filename = "../submission/pos-neg-patterns_PerMod.png", width = 6, height = 6, units = "in", res = 500)
plot_expression_profiles(UHRBHR_expr_filtered, modules$modules)
# dev.off()


#### Functional enrichment ####
# Q4. Can you explain the enrichment result shown on the resulting plot

enrichment <- bio_enrich(modules$modules, organism="hsapiens")
# if not pruning the version number --> Error in bio_enrich(modules$modules, organism = "hsapiens") : 
#                                           could not find function "bio_enrich"

# Module 1
enrichment$result %>% dplyr::filter(query == 1)
# Module 2
enrichment$result %>% dplyr::filter(query == 2)

### Plot 4
p = plot_enrichment(enrichment)
p
# htmlwidgets::saveWidget(p, file = "../submission/enrichment_PerMod.html", selfcontained = TRUE)

