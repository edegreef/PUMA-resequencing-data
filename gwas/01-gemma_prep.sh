#!/bin/bash

#Getting the files ready for running Purple Martin GWAS with Gemma
input_vcf=/scratch/user/edegreef/reseq/filter_scaffolds/imputed.qualityfiltered.geobirds89.recode.vcf #input file
phenotype=spring_depart #input file
output_plink=SD.allchrom.89
type=SD_allchrom_89

#converting to plink format and adding phenotype
/general/home/edegreef/plink/plink --vcf $input_vcf --make-bed --allow-extra-chr --out $output_plink
awk '{print $1, $2, $3, $4, $5}' $output_plink.fam > temp.fam
paste temp.fam $phenotype > $output_plink.fam
rm temp.fam

mkdir $type
mv $output_plink.bed $type/$output_plink.bed
mv $output_plink.bed $type/$output_plink.bim
mv $output_plink.bed $type/$output_plink.fam

