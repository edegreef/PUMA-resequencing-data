#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=gem_prep
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 


#indexing genome
/scratch/user/edegreef/gemtools-1.7.1-i3/gemtools index -i ref_genome/PUMA2_03022020.fasta

#creating mappability file 
/scratch/user/edegreef/gemtools-1.7.1-i3/gem-mappability -I PUMA2_03022020.gem -l 150 -o PUMA2_03022020_mappability_150