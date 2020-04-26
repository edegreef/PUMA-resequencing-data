#!/bin/bash

#Filtering Purple Martin SNPs round 2: making multiple files containg autosomes, population birds, & geolocator birds

#input files:
input_vcf=imputed.qualityfiltered
sex_scaf=fst_sex_linked_chrom_pos
pop_indiv=list_4pop_indiv
geo_birds=list_geobirds_indiv

#Creating vcf of autosomes (excluding Z & W chromosomes)
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
vcftools --vcf $input_vcf.recode.vcf --exclude-positions $sex_scaf --recode --recode-INFO-all --out $input_vcf.Zout

#Subsetting vcf for pop birds and geolocator birds
#including all chromosomes
vcftools --vcf $input_vcf.recode.vcf --keep $pop_indiv --recode --recode-INFO-all --out $input_vcf.pop
vcftools --vcf $input_vcf.recode.vcf --keep $geo_indiv --recode --recode-INFO-all --out $input_vcf.geobirds

#including only autosomes
vcftools --vcf $input_vcf.ZWout.recode.vcf --keep $pop_indiv --recode --recode-INFO-all --out $input_vcf.ZWout.pop
vcftools --vcf $input_vcf.ZWout.recode.vcf --keep $geo_indiv --recode --recode-INFO-all --out $input_vcf.ZWout.geobirds
