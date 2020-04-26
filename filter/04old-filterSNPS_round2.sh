#!/bin/bash

#Filtering Purple Martin SNPs round 2: making multiple files containg autosomes, Z chrom, W chrom, population birds, & geolocator birds

#input files:
input_vcf=imputed.qualityfiltered
Z_scaf=Z_scaffolds_MF0.04_chrom_pos
W_scaf=W_scaffolds_M0.7_F0.8_chrom_pos
pop_indiv=list_4pop_indiv
geo_indiv=list_geobirds_indiv

#Step1: Creating vcf of autosomes (excluding Z & W chromosomes)
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
vcftools --vcf $input_vcf.recode.vcf --exclude-positions $Z_scaf --recode --recode-INFO-all --out $input_vcf.Zout
vcftools --vcf $input_vcf.Zout.recode.vcf --exclude-positions $W_scaf --recode --recode-INFO-all --out $input_vcf.ZWout
rm $input_vcf.Zout.recode.vcf

#Step2: Extracting just the Z chromosome and the W chromosome (separately)
vcftools --vcf $input_vcf.recode.vcf --positions $Z_scaf --recode --recode-INFO-all --out $input_vcf.Zonly
vcftools --vcf $input_vcf.recode.vcf --positions $W_scaf --recode --recode-INFO-all --out $input_vcf.Wonly

#Step3: Subsetting vcf for pop birds (N=71) and geolocator birds (N=93)
#including all chromosomes
vcftools --vcf $input_vcf.recode.vcf --keep $pop_indiv --recode --recode-INFO-all --out $input_vcf.pop
vcftools --vcf $input_vcf.recode.vcf --keep $geo_indiv --recode --recode-INFO-all --out $input_vcf.geobirds
#including only autosomes
vcftools --vcf $input_vcf.ZWout.recode.vcf --keep $pop_indiv --recode --recode-INFO-all --out $input_vcf.ZWout.pop
vcftools --vcf $input_vcf.ZWout.recode.vcf --keep $geo_indiv --recode --recode-INFO-all --out $input_vcf.ZWout.geobirds
#including only Z chromosome
vcftools --vcf $input_vcf.Zonly.recode.vcf --keep $pop_indiv --recode --recode-INFO-all --out $input_vcf.Zonly.pop
vcftools --vcf $input_vcf.Zonly.recode.vcf --keep $geo_indiv --recode --recode-INFO-all --out $input_vcf.Zonly.geobirds
