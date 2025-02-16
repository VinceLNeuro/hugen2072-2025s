Step 4. 
* IBD KIND: 
```
Excluding 1,339 SNPs (non-autosomes or non-selection)
Excluding 35 SNPs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
    # of samples: 77
    # of SNPs: 1,926
```

=== Get genetic relatedness ===

To calculate expected relatedness:
* Remove duplicates from the pedigree (scanAnnot > family info): 43 uniq left
* Add (known ID but not included in the pedigree) parent info (n=1, in family 1341) to the pedigree data
* Fix the half-parental missing in family 58 (add n=1)
* Subfamily structure: 1341 and 1362 have been split into 2 trios (the new families are 1341-2 and 1362-2)

=== calculate expected relatedness ===

* no inbreeding parents
* expected/reported relationship table
    - PO  U
    - 30 15
    
=== Plotting ===

* Identify unmatched relationship-pairs
* **Fixed ScanID 355 and 356 (who are actually unrelated => give one of them a new family ID)**
* **Fixed Sample-swap (PO-U mismatch) of 200016815 (296,297) & 200071490 (298,299)
    - 2 subjects, but 4 scans
    
---

Step 5. PCA

1. Exclude SNPs from a few problematic loci: 2q21 (LCT), HLA, 8p23, and 17q21.31: no snps removed
2. Exclude duplicates for a subject (only keep the lowest missing.e1 one): kept 43 subjects 
3. LD pruning -> 177 independent SNPs
4. conduct PCA
5. visualization

---

Step 6. Duplicate Sample Discordance to Assess *Probe Quality*

1. (using missing.n1 as threshold) 25 SNPs has at least 1 discordant pairs
2. (using missing.n2 as threshold) 21 SNPs were filtered (≥ 1 discordant pairs)

---

Step 7. HWE

1. Exclude Scans: populations other than CEU, nonfounders, and duplicates (exclude_n = 60)

    - 77 scans total; excluding 60 scans => **17 kept**
    - matches sample sizes for autosomes

2. After removing the monomorphic & missing SNPs, 1785 autosome SNPs left and 700 sex chromosome SNPs left (in total 2485 SNPs)
    


