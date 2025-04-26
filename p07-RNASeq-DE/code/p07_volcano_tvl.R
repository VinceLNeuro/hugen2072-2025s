# call in RStudio, due to .Rproj file
# setwd("p07-RNASeq-DE/code/")
rm(list=ls())
library(EnhancedVolcano)
library(tidyverse)
library(latex2exp)

# # Get the summary table for **all genes** with FDR values
# tab_summary_allGenes = 
#     topTags(et12,    # exactTest result
#     n = Inf,         # all genes
#     adjust.method = "BH", # FDR correction
#     p.value = 1)     # all genes
# 
# # Checks
# tab_summary_allGenes %>% nrow() # "654 genes are left"
# tab_summary_allGenes %>% as.data.frame() %>% filter(FDR<0.05) %>% nrow() # 77 DEGs
# 
# # Store the summary table
# df_summary_allGenes = rownames_to_column(as.data.frame(tab_summary_allGenes), var = "ENSG")
# nrow(df_summary_allGenes)
# write_csv(df_summary_allGenes, file="../output/3_DiffGeneExpr_result/table_AllGenes_fdr.csv")


# Load data
summayTab = read_csv("../output/3_DiffGeneExpr_result/table_AllGenes_fdr.csv", col_names = TRUE) %>% as.data.frame()
rownames(summayTab) = summayTab$ENSG
summayTab = summayTab %>% arrange(FDR) %>% select(-ENSG) #sorted by FDR
head(summayTab)

# Create lookup table for plotting
keyval = rep('grey', nrow(summayTab))      #value=col
names(keyval) = rep('NS', nrow(summayTab)) #names=annot
#sig-up
keyval[which(summayTab$FDR < 0.05 & summayTab$logFC > 0)] = 'red'
names(keyval)[which(summayTab$FDR < 0.05 & summayTab$logFC > 0)] = 'Upregulated'
#sig-down
keyval[which(summayTab$FDR < 0.05 & summayTab$logFC < 0)] = 'blue'
names(keyval)[which(summayTab$FDR < 0.05 & summayTab$logFC < 0)] = 'Downregulated'
#check the status of keyval
unique(names(keyval))
unique(keyval)
keyval[1:20]

# Create the volcano plot
topN = 15
vol = EnhancedVolcano(summayTab,
                      lab = rownames(summayTab), #all labels
                      x = "logFC",
                      y = "FDR", ylab = TeX("$-Log_{10}FDR$"),
                      xlim = c(-12,12),
                      ylim = c(0,6),
                      title = 'DE mRNAs comparing HBR and UHR',
                               # UHR = cancer cell lines
                               # HBR = brains of 23 Caucasians of varying age but mostly 60-80 years old
                      subtitle = paste0("Top ",topN," labeled"),
                      caption = paste0('Total = ',nrow(summayTab),' ENSG\n
                                        Total DE = ',nrow(summayTab %>% filter(FDR<0.05)),' ENSG'),
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      cutoffLineCol = 'black', cutoffLineType = 'dashed',
                      selectLab = rownames(summayTab)[1:topN],
                      colCustom = keyval,
                      
                      drawConnectors = TRUE, max.overlaps = Inf, maxoverlapsConnectors = Inf, min.segment.length = 0.2, # minimum length without an arrow
                      pointSize = 2, legendIconSize = 2,
                      labSize = 2,
                      axisLabSize = 12, titleLabSize = 15, legendLabSize = 12, subtitleLabSize = 12, captionLabSize = 12,
                      legendPosition = 'bottom',
                      border = 'full', borderWidth = 0.5)
# vol

# Save the volcano plot
# ggsave(path = '../output/3_DiffGeneExpr_result/', 'vol_HBR-UHR.png',
#        plot = vol,
#        device = 'png', dpi = 500,
#        width = 8, #if full size use 12,8
#        height = 8)
