#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=msmc_boot_run_pmref
#SBATCH --time=12:00:00
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
boot_list=PMREF_bootlist

while read boot
do
	/scratch/user/edegreef/msmc/msmc2 -t 6 -o boots/output/PMREF.min1MB.$boot.msmc2 boots/$boot/bootstrap_multihetsep*
done < $boot_list


