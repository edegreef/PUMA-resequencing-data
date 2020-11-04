#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=msmc_run
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=8G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 


#use 32G for larger job
cd /scratch/user/edegreef/msmc

#running msmc with all multihetsep files for  one individual
/scratch/user/edegreef/msmc/msmc2 -t 6 -o msmc_output/PM150.min1MB.msmc2 multihetsep_PM150_min1MB/PM*


