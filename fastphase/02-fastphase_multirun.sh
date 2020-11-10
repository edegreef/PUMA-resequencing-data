#!/bin/sh
## author: kira delmore, modified by evelien de greef
## date: mar 2018 - modified nov 2020

#need to have the fastphase_subpoplabels.inp file in the working directory (space delimited file with pop info)

list=list_scaf_all
job="fastphase"

mkdir fastphase_input_all
mkdir fastphase_output_all

while read scaffold
do

echo "#BSUB -L /bin/bash              
#BSUB -J fastphase."$scaffold"
#BSUB -n 8
#BSUB -R "span[ptile=8]" 
#BSUB -R "rusage[mem=2500]"
#BSUB -M 2500
#BSUB -W 24:00
#BSUB -o stdout.%J."$scaffold".fastphase.err
#BSUB -e stderr.%J."$scaffold".fastphase.out

cd /scratch/user/edegreef/reseq/fastphase

#convert vcf to fastphase input formats (geno and markers file)
perl vcf2fastPHASE.pl split_vcf_all/imputed.qualityfiltered.87."$scaffold".recode.vcf fastphase_input_all/"$scaffold".87.geno fastphase_input/"$scaffold".87.markers 87

#run fastphase using the geno files
/scratch/user/edegreef/reseq/fastphase/fastPHASE -u -C25 -H100 -K10 -ofastphase_output_all/"$scaffold".87 fastphase_input_all/"$scaffold".87.geno


" > $job.$scaffold.lsf

bsub < $job.$scaffold.lsf

done < $list

mkdir job_logs
mv fastphase.scaffold* job_logs
mv *fastphase.out job_logs
mv *fastphase.err job_logs
