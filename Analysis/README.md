# Downstream analysis
We create a vcf file containing small and structral variants annotated with their AF and predicted functional impact. Here we are going to go through how to explore those variants.


If we one want to work on the SV, we can filter the vcf using this command from gatk
```
gatk SelectVariants -V <input.vcf> --select-type-to-include INDEL --min-indel-size 50 -O output.vcf
```

# Visualize SVs
We benchmarked the 
###samplot

## plot a SV
Lets first get high impact variant
```
zgrep "|HIGH|" results/LR_calling/variants/GG/cattle_taurus_10.ERR7091271.chr25.ont.minimap2/annotated/merged.vcf.gz 
```

lets visualize 
```
snakemake -np ../results/LR_calling/samplot/DEL_25_41669173_41669210.png
```


### ribbon

## Population Frequency
The Great Genotyper estimate the genotypes for the found variants in 464 cattle samples. We used bcftools to calculate some metrics using the genotypes. Variants at [path] are tagged by the following metadata. 
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

To plot population histogram for specific SV
```
 grep <VID> | calculate histogram
```

To calculate AF for specific breeds
```
```

To view sample names having a specific variants
```
```




## Variant Effect predictor 
We used VEP to predict the effect of the variants, the output annotated vcf is at [path] and the reprot is at [path]. The following figgure summarizes the annotations produced by VEP. for more information visit their website
![VEP](https://uswest.ensembl.org/info/genome/variation/prediction/consequences.jpg)

To view novel High impact variants 
```
bcftools view  -q 0.001 <variant.vcf> | grep "|HIGH|"
```

To filter variants affecting specific gene
```
bcftools view  -q 0.001 <variant.vcf> | grep "|<GENENAME>|"
```
