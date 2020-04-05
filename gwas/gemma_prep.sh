#!/bin/bash

#Getting the input files ready for running Purple Martin GWAS with Gemma

input_vcf=/scratch/user/edegreef/reseq/filter_scaffolds/imputed.qualityfiltered.geobirds.recode.vcf
output_plink=imputed.qualityfiltered.geobirds
type=allchrom

#input lists of phenotype
SA_mig=spring_arrive
SD_mig=spring_depart
FA_mig=fall_arrive
FD_mig=fall_depart

#converting to plink format
/general/home/edegreef/plink/plink --vcf $input_vcf --make-bed --allow-extra-chr --out $output_plink
awk '{print $1, $2, $3, $4, $5}' $output_plink.fam > temp.fam

#setting up spring arrival files
mkdir SA_$type
cp $output_plink.bed SA_$type/SA.$type.bed
cp $output_plink.bim SA_$type/SA.$type.bim
paste temp.fam $SA_mig > SA_$type/SA.$type.fam

#setting up spring depart files
mkdir SD_$type
cp $output_plink.bed SD_$type/SD.$type.bed
cp $output_plink.bim SD_$type/SD.$type.bim
paste temp.fam $SD_mig > SD_$type/SD.$type.fam

#setting up fall arrival files
mkdir FA_$type
cp $output_plink.bed FA_$type/FA.$type.bed
cp $output_plink.bim FA_$type/FA.$type.bim
paste temp.fam $FA_mig > FA_$type/FA.$type.fam

#setting up fall depart files
mkdir FD_$type
cp $output_plink.bed FD_$type/FD.$type.bed
cp $output_plink.bim FD_$type/FD.$type.bim
paste temp.fam $FD_mig > FD_$type/FD.$type.fam

rm temp.fam