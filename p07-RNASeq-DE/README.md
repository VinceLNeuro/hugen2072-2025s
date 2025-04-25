# Project 07 - RNA-seq for DE

__Author__: Tianze (Vincent) Luo \
__Date__: 2025-04-24


## 0. Background

* RNA-seq data on a single chromosome in multiple experiments

* Files: `/ix1/hugen2072-2025s/p7/`

* Interactive session command: `crc-interactive --teach -a hugen2072-2025s -t 4:00:00`

<br>

## 1. `cutadapt`: trim adapter sequence

<br>

## 2. `STAR`: Aligning the reads 

- Perform downstream analysis using `*ReadsPerGene.out.tab`

```sh
[til177@teach-cpu-n0 2_STAR]$ ll *ReadsPerGene.out.tab
# -rw-r--r-- 1 til177 mmarazita 1.5M Apr 24 23:08 HBR_Rep1ReadsPerGene.out.tab
# -rw-r--r-- 1 til177 mmarazita 1.5M Apr 24 23:08 HBR_Rep2ReadsPerGene.out.tab
# -rw-r--r-- 1 til177 mmarazita 1.5M Apr 24 23:08 UHR_Rep1ReadsPerGene.out.tab
# -rw-r--r-- 1 til177 mmarazita 1.5M Apr 24 23:08 UHR_Rep2ReadsPerGene.out.tab
```

<br>

## 3. DGE analysis in R

```sh
module load gcc/12.2.0 r/4.4.0
Rscript --vanilla p07_DE_tvl.R > ../output/DiffGeneExpr.log
```

- see `output/DiffGeneExpr.log` and `code/Rplots.pdf` for details.
