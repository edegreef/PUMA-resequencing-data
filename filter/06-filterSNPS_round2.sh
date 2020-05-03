#!/bin/bash

#Filtering Purple Martin SNPs round 2: making multiple files containg autosomes, population birds, & geolocator birds

#input files:
input_vcf=imputed.qualityfiltered.LD
sex_scaf=fst_sex_linked_chrom_pos_LD
pop_indiv=list_4pop_indiv
geo_indiv=list_geobirds_indiv

#Creating vcf of autosomes (excluding Z & W chromosomes)
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
vcftools --vcf $input_vcf.recode.vcf --exclude-positions $sex_scaf --recode --recode-INFO-all --out $input_vcf.ZWout

#Subsetting vcf for pop birds and geolocator birds
#including all chromosomes
vcftools --vcf $input_vcf.recode.vcf --keep $pop_indiv --recode --recode-INFO-all --out $input_vcf.pop
vcftools --vcf $input_vcf.recode.vcf --keep $geo_indiv --recode --recode-INFO-all --out $input_vcf.geobirds

#including only autosomes
vcftools --vcf $input_vcf.ZWout.recode.vcf --keep $pop_indiv --recode --recode-INFO-all --out $input_vcf.ZWout.pop
vcftools --vcf $input_vcf.ZWout.recode.vcf --keep $geo_indiv --recode --recode-INFO-all --out $input_vcf.ZWout.geobirds


###adding on section for files that are excluding the 4 sample outliers (PM65, 139, 174b, 182)
vcftools --vcf $input_vcf.recode.vcf --keep list_4pop_indiv_67 --recode --recode-INFO-all --out $input_vcf.pop67
vcftools --vcf $input_vcf.recode.vcf --keep list_geobirds_indiv_89 --recode --recode-INFO-all --out $input_vcf.geobirds89

#including only autosomes - geo
vcftools --vcf $input_vcf.ZWout.recode.vcf --keep list_4pop_indiv_67 --recode --recode-INFO-all --out $input_vcf.ZWout.pop67
vcftools --vcf $input_vcf.ZWout.recode.vcf --keep list_geobirds_indiv_89 --recode --recode-INFO-all --out $input_vcf.ZWout.geobirds89