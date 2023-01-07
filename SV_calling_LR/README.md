# Intro
We are going to learn how to study structral variants(SV) in cattle genomes using long reads from Pacbio and Oxford nanpore. First, We are going to discover and phase  small variants, and SVs using state of the art tools. After that, We are going to calculate population allele frequencies using novel population genotyper tool(The Great Genotyper).  Lastly, We are going to functionally annotate the variants to facilitatle studying the functional impact of the SVs.

We chose sequencing datasets from  haplotype-resolved assembly project(PRJEB42335) of Nellore and Brown_Swiss cross for two reasons: First, We can create a gold standard benhcmark by calling the variants from the haplotype-resolved assemblies which is considered to be the accurate method Figure 1. Second, the sample was heavily sequenced using illumina, pacbio(HIFI), and oxford nanpore which allows us to compare the results of differnet metohds.  Lastly, we  created a downsampled the data for the sake of the workshope. We are going to focus on chromsome 25 only, and we are going to calculate the AF in 30 samples(scalable to 4000 samples). 


|![compare](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-019-1828-7/MediaObjects/13059_2019_1828_Fig2_HTML.png?as=webp)|
|:--:|
|Figure 1: Comparison of different methods
Ref:  Mahmoud M, Gobet N, Cruz-Dávalos DI, Mounier N, Dessimoz C, Sedlazeck FJ. Structural variant calling: The long and the short of it. Genome Biology. 2019 Nov 20;20(1):246. 
|



# Data Description 
Test Input data for the workshop can be downloaded from the following link. 
Data description is
| file        |  Description  |
|:-------------:|:-------------|
| ARS-UCD1.2_Btau5.0.1Y.25.fa | Chromosome 25 from ARS-UCD1.2 genome|
| ARS-UCD1.2_Btau5.0.1Y.25.rmsk.bed.gz | Repeats annotation of chromosome 25|
| ARS-UCD1.2_Btau5.0.1Y.25.gff.gz | genes annotation of chromosome 25|
| goldstandard/callset_filered.25.vcf.gz | gold standard for variant calling created using the haployte resolved assemblies|
| goldstandard/callset_filered.25.bed.gz | the gold standard in bed format for plotting|
| ERR5043144.chr25 | Hifi reads from sample SAMEA10017982 that maps to chromsome 25|
| ERR7091271.25.fastq.gz | ONT reads from sample SAMEA10017982 that maps to chromsome 25|
| cattle_taurus_10 |  Folder contains Kmer indexes of 10  taurus samples|
| cattle_taurus_10/graph.desc.tsv |  file contains the Biosample ids|
| cattle_indicus_10 |  Folder contains Kmer indexes of 10  indicuis samples|
| cattle_indicus_10/graph.desc.tsv |  file contains the Biosample ids|
| cattle_bostgroup_10 |  Folder contains Kmer indexes of 10  boison samples|
| cattle_bosgroup_10/graph.desc.tsv |  file contains the Biosample ids|


# Workflow:

We created a snakemake script to wrap all the commands in the tutorial. The workflow(summarized in Figure 2) has the following steps: 
  1. map using minmap2
  2. call small variants using clair3
  3. Phase small variants using longshot
  4. split the reads into two haplotypes
  5. call SV using pbsv, sniffles, cuteSV
  6. merge the small and structrual variants
  7. annotate the vcf using Variant effect predictor
  8. calculate AF using the great genotyper

![worflow_dag](dag.png)



 
# How to Run
1. add your samples names and sequencing technology to sample_table.csv
2. add samples files to subsample_table.csv
3. add outputfolder to config.yaml
4. run using snakemake -j16 --use-conda
