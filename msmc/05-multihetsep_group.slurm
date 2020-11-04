#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=msmc_multihetsep_group5
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

cd /scratch/user/edegreef/msmc

#code for generating multihetsep outputs for group of 5 samples

list=list_scaf_min1MB

sample1=PM190_AB_min1MB
sample2=PM201_AB_min1MB
sample3=PM196_AB_min1MB
sample4=PM135_AB_min1MB
sample5=PM105_AB_min1MB


mkdir multihetsep_AB_select5

module load Python/3.8.2-GCCcore-9.3.0

while read scaffold
do
	msmc-tools/generate_multihetsep.py --mask mask_output_$sample1/$sample1.mask.$scaffold.bed.gz --mask mask_output_$sample2/$sample2.mask.$scaffold.bed.gz --mask mask_output_$sample3/$sample3.mask.$scaffold.bed.gz --mask mask_output_$sample4/$sample4.mask.$scaffold.bed.gz --mask mask_output_$sample5/$sample5.mask.$scaffold.bed.gz --mask /scratch/user/edegreef/gem_prep/PUMA2_03022020_mappability_150.bed mask_output_$sample1/$sample1.$scaffold.vcf.gz mask_output_$sample2/$sample2.$scaffold.vcf.gz mask_output_$sample3/$sample3.$scaffold.vcf.gz mask_output_$sample4/$sample4.$scaffold.vcf.gz mask_output_$sample5/$sample5.$scaffold.vcf.gz > multihetsep_AB_select5/$scaffold.multihetsep.txt
done < $list