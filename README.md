# workshop_12Jan_2023

# 1) Intro

We are going to learn how to study structural variants (SV) in cattle
genomes using long reads from Pacbio and Oxford Nanopore.  We are
going to discover and phase small variants and SVs using state of the
art tools. After that, we are going to calculate population allele
frequencies using a novel population genotyper tool (The Great
Genotyper).  Lastly, we are going to functionally annotate the
variants to facilitate studying the functional impact of the SVs.

## Aims

1. learn how to call and phase SVs using the different SV callers
2. compare the performance of different methods
3. learn how to calculate population allele frequency for cattle SVs
4. learn how to predict the impact of SV on genes

# 2) Significance 

* Novel way of calculating Population AF without needing databases like genomAD.
* Variants will be richly annotated with AF and predicted functional impact making it perfect for studying functional impact of the SVs.
* Workflow is benchmarked on cattle data and achieving accuracy of ~ 90%.
* Workflow is implemented using snakemake to make it simpler to run it afterward with your data and tweak it as you want.

# Tutorial is in SV_calling_LR
