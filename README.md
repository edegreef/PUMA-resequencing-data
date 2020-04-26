# PUMA-SNPs

This is a repository for scripts I used in analyzing Purple Martin (*Progne subis*) resequencing data.
These scripts were executed using the Texas A&M High Preformance Research Computing resources (https://hprc.tamu.edu/).

# Filtering SNPs 
[:file_folder:](https://github.com/edegreef/PUMA-SNPs/tree/master/filter)

Primary input file:
* SNPs in vcf format

Separate files:
* scaffolds_read_depths.csv (read depth information of the samples along the reference genome)
* filter.hets.R (included in filter folder)
* run_script.lsf (included in filter folder)

1. Create list of CHROM and POS info from each SNP (used in downstream filtering steps for particular scaffolds)
2. ID scaffolds that are "weird" (with abnormally high read depths), Z-linked, and W-linked, using the list_all from step1 and scaffolds_read_depths file. Lists created from this step are used in filtering steps (step3 & 4)
UPDATE::: ended up using FST for male/females to filter out sex-linked scaffolds.
3. Filter SNPs based on quality (indels, minimum allele frequency, missing genotypes, Hardy-Weinberg Equilibrium, triallelic, and weird/high read depths)
4. Filter SNPs to create files for autosomes, Z chromosome, W chromosome, and separate files for geolocator birds and population birds.



# FST
in prog

# PCA
in prog

# GWAS
in prog
