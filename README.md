# gwas
Genome-wide association study

A small genome-wide association study on simulated data.\
The file gwas_population.vcf contains the genotypes of 1000 individuals at 10,000 SNPs.\
The file gwas_phenotypes.txt describes which of those 1000 individuals have disease D (a value of 1 
indicates that the individual has the disease).

This is a program that reads these two files and computes, for each SNP, the 
p-value of the association between the SNP’s genotype and the disease D, using a 
standard chi-squared test.
