# Project 06: Structual Variant Calling Pipeline Report

__Author__: Tianze (Vincent) Luo <br>
__Date__  : 2025-04-21 <br>

__Table of Contents__
- [Project 06: Structual Variant Calling Pipeline Report](#project-06-structual-variant-calling-pipeline-report)
  - [0. Background](#0-background)
  - [1. Extract chr22 alignments (`.BAM`) \&\& index](#1-extract-chr22-alignments-bam--index)
  - [2. Manta caller (used single tool in this pipeline)](#2-manta-caller-used-single-tool-in-this-pipeline)
  - [3. Convert chr22 alignments `.bam` to indexed `.cram`](#3-convert-chr22-alignments-bam-to-indexed-cram)
  - [4. Investigate with IGV](#4-investigate-with-igv)


<br>

## 0. Background

- Interactive session command: `srun -N 1 --pty bash`

<br>

- Start `.cram` file: `p6/NA12778.final.cram`
  
  - Alternative format for `.bam`
  
  - __NA12778__: Female resident of Utah with Northern and Western European–associated ancestry [1000 Genome Project]

- Reference sequence: `p6/GRCh38_full_analysis_set_plus_decoy_hla.fa`

<br>


## 1. Extract chr22 alignments (`.BAM`) && index

Will call SVs on chromosome 22 only


## 2. Manta caller (used single tool in this pipeline)

manta SV discovery caller based on __discordant pairs (PE)__ and __split reads (SR)__.

`manta` workflow:
- `configManta.py`
- `manta_test/runWorkflow.py`
- variant `.VCF` -> `.BED`
  - Extract <mark>DEL/DUP > 1 kb<mark> (`.bed`)


## 3. Convert chr22 alignments `.bam` to indexed `.cram`

<br>

## 4. Investigate with IGV