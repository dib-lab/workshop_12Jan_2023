# Analysis
Snakemake workflow will create vcf files containing small and structral variants annotated with their AF and predicted functional impact. Here we are going to go through how to navigate those vcfs.

## General analysis
Calculate the total number of variants
```
bcftools view <input.vcf> |grep -vP '^#' |wc -l
```

Calcualte number of variants per chromsome
```
bcftools view <input.vcf> |grep -vP '^#' |cut -f1 |uniq -c
```

Filter SV
```
gatk SelectVariants -V <input.vcf> --select-type-to-include INDEL --min-indel-size 50 -O output.vcf
```

## Visualize SVs



## Population Frequency
The Great Genotyper estimate the genotypes for the found variants in 500 cattle samples. We used bcftools to calculate some metrics using the genotypes. Variants at [path] are tagged by the following metadata. 
| Metadata      | Description |
| -- |:-----------:|
| AC | Allele count in genotypes|
| AC_Het | Allele counts in homozygous genotypes|
| AC_Hom | Allele counts in heterozygous genotypes|
| AC_Hemi | Allele counts in hemizygous genotypes|
| AF | Allele frequency |
| MAF | Minor Allele frequency |
| NS | Number of samples with data   |
| AN | Total number of alleles in called genotypes |
| HWE | Hardy-Weinberg equilibrium |
| ExcHet | Test excess heterozygosity; 1=good, 0=bad |


Bcftools is very helpful in filtering vcf files using the variants metadata. For example, We can query the novel varaints using the following command
```
bcftools view  -Q 0.001 <variant.vcf> -o novel.vcf 
```
or 
```
bcftools view  -i "INFO/AF <= 0.0001"  <variant.vcf> -o novel.vcf 
```




## Variant Effect predictor 
We used VEP to predict the effect of the variants, the output annotated vcf is at [path] and the reprot is at [path]. The following figgure summarizes the annotations produced by VEP. for more information visit their website
![VEP](https://uswest.ensembl.org/info/genome/variation/prediction/consequences.jpg)

To view novel High impact variants 
```
bcftools view  -q 0.001 <variant.vcf> | grep "|HIGH|"
```

