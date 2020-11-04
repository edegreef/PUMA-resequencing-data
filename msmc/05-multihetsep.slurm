#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=multihetsep
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

cd /scratch/user/edegreef/msmc

#code for generating multihetsep outputs for individal sample
sample1=PM150_min1MB
list=list_scaf_min1MB

mkdir multihetsep_$sample1

module load Python/3.8.2-GCCcore-9.3.0

while read scaffold
do
	msmc-tools/generate_multihetsep.py --mask mask_output_$sample1/$sample1.mask.$scaffold.bed.gz --mask /scratch/user/edegreef/gem_prep/PUMA2_03022020_mappability_150.bed mask_output_$sample1/$sample1.$scaffold.vcf.gz > multihetsep_$sample1/$sample1.$scaffold.multihetsep.txt
done < $list
echo "PM150 run done"
