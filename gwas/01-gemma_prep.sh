#!/bin/bash

#Getting the input files ready for running Purple Martin GWAS with Gemma

input_vcf=/scratch/user/edegreef/reseq/filter_scaffolds/imputed.qualityfiltered.geobirds.recode.vcf
file=SA_allchrom
migdata=spring_arrive

#converting to plink format
/general/home/edegreef/plink/plink --vcf $input_vcf --make-bed --allow-extra-chr --out $file

awk '{print $1, $2, $3, $4, $5}' $file.fam > temp.fam
paste temp.fam $migdata > $file.2.fam 
rm $file.fam temp.fam
mv $file.2.fam $file.fam

mkdir temp
mv $file* ./temp
mv temp $file