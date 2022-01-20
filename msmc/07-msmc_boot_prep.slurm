#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=msmc_boot
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=8G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 


cd /scratch/user/edegreef/msmc

module load Python/3.8.2-GCCcore-9.3.0

msmc-tools/multihetsep_bootstrap.py boot_PM141 multihetsep_PM141_min1MB/PM141_min1MB.*.multihetsep.txt
msmc-tools/multihetsep_bootstrap.py boot_PM150 multihetsep_PM150_min1MB/PM150_min1MB.*.multihetsep.txt
msmc-tools/multihetsep_bootstrap.py boot_PM223 multihetsep_PM223_min1MB/PM223_min1MB.*.multihetsep.txt
msmc-tools/multihetsep_bootstrap.py boot_PMREF multihetsep_PMREF_min1MB/PMREF_min1MB.*.multihetsep.txt

