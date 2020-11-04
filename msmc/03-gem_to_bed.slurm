#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=gem_to_bed
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

cd /scratch/user/edegreef/gem_prep

/scratch/user/edegreef/gemtools-1.7.1-i3/gem-2-wig -I PUMA2_03022020.gem -i PUMA2_03022020_mappability_150.mappability -o PUMA2_03022020_mappability_150

/scratch/user/edegreef/tools/wigToBigWig PUMA2_03022020_mappability_150.wig PUMA2_03022020_mappability_150.sizes PUMA2_03022020_mappability_150.bw

/scratch/user/edegreef/tools/bigWigToBedGraph PUMA2_03022020_mappability_150.bw  PUMA2_03022020_mappability_150.bedGraph

/scratch/user/edegreef/tools/bedGraphTobed PUMA2_03022020_mappability_150.bedGraph PUMA2_03022020_mappability_150.bed 0.3