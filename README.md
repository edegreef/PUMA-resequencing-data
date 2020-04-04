# PUMA-SNPs

This is a repository for scripts I used in analyzing Purple Martin (*Progne subis*) resequencing data.
These scripts were executed using the Texas A&M High Preformance Research Computing resources (https://hprc.tamu.edu/).

# Filtering SNPs 
[:file_folder:](https://github.com/edegreef/PUMA-SNPs/tree/master/filter)

Primary input file is SNPs in vcf format. I also had scaffolds_read_depths (for X reason), and a few other files for running R

1. Create list of CHROM and POS info from each SNP (used in downstream filtering steps for particular scaffolds)
2. Use list from step1, and info from scaffolds_read_depths file (separate) to ID weird scaffolds with higher read depths to filter
3. Filter SNPs based on quality (indels, minimum allele frequency, missing genotypes, Hardy-Weinberg Equilibrium, triallelic, and weird/high read depths)
4. Identify sex-linked scaffolds
5. Filter SNPs to create files for autosomes, Z chromosome, W chromosome, and separate files for geolocator birds and population birds.
