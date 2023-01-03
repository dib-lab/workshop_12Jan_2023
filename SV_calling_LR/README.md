#Summary
Workflow to call, phase and annotate small, and SV from long reads. AF for SV will be annotated using novel population genotyping technique agains 500 cattle sample. 
Workflow:
  1. map using minmap2
  2. call small variants using clair3
  3. Phase small variants using longshot
  4. split the reads into two haplotypes
  5. call SV using pbsv, sniffles, cuteSV
  6. merge the small and structrual variants
  7. annotate the vcf using Variant effect predictor
  8. calculate AF using the great genotyper
#How to Run
1. add your samples names and sequencing technology to sample_table.csv
2. add samples files to subsample_table.csv
3. add outputfolder to config.yaml
4. run using snakemake -j16 --use-conda
