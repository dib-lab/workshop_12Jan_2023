# Summary

Workflow to call, phase and annotate small, and SV from long reads. AF for SV will be annotated using novel population genotyping technique agains 464 cattle samples. 

Workflow:
  1. map using minmap2
  2. call small variants using clair3
  3. Phase small variants using longshot
  4. split the reads into two haplotypes
  5. call SV using pbsv, sniffles, cuteSV
  6. merge the small and structrual variants
  7. annotate the vcf using Variant effect predictor
  8. calculate AF using the great genotyper

![worflow_dag](dag.pdf)

# Data needed
Test Input data for the workshop can be downloaded from the following link. 
Data description is
| file        |  Description          |
|:-------------:|:-------------|
| ARS-UCD1.2_Btau5.0.1Y.25.fa | Chromosome 25 from ARS-UCD1.2 genome|
| ERR5043144.chr25 | Hifi reads from sample SAMEA10017982 that maps to chromsome 25|
| ERR7091271.25.fastq.gz | ONT reads from sample SAMEA10017982 that maps to chromsome 25|
| population_index_1 |  Folder contains Kmer indexes of 225 samples|
| population_index_1/graph.desc.tsv |  file contains 225 Biosample ids|
| population_index_2 |  Folder contains Kmer indexes of 203 samples|
| population_index_2/graph.desc.tsv |  file contains 203 Biosample ids|


 
# How to Run
1. add your samples names and sequencing technology to sample_table.csv
2. add samples files to subsample_table.csv
3. add outputfolder to config.yaml
4. run using snakemake -j16 --use-conda
