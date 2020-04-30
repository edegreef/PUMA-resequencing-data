#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=gemma_run_bslmm
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

cd /scratch/user/edegreef/gemma/SD_allchrom_89
file=SD.allchrom.89

#making relatedness matrix
/general/home/edegreef/gemma-0.98.1-linux-static -bfile $file -gk 1 -miss 1 -maf 0 -notsnp -o relatedness_matrix

#running BSLMMs
/general/home/edegreef/gemma-0.98.1-linux-static -bfile $file -k output/relatedness_matrix.cXX.txt -miss 1 -maf 0 -notsnp -bslmm 1 -w 5000000 -s 20000000 -rpace 100 -wpace 1000 -hmin 0 -hmax 1 -rmin 0 -rmax 1 -o $file.bslmm_1
echo "done first gemma run"

/general/home/edegreef/gemma-0.98.1-linux-static -bfile $file -k output/relatedness_matrix.cXX.txt -miss 1 -maf 0 -notsnp -bslmm 1 -w 5000000 -s 20000000 -rpace 100 -wpace 1000 -hmin 0 -hmax 1 -rmin 0 -rmax 1 -o $file.bslmm_2
echo "done second gemma run"

/general/home/edegreef/gemma-0.98.1-linux-static -bfile $file -k output/relatedness_matrix.cXX.txt -miss 1 -maf 0 -notsnp -bslmm 1 -w 5000000 -s 20000000 -rpace 100 -wpace 1000 -hmin 0 -hmax 1 -rmin 0 -rmax 1 -o $file.bslmm_3
echo "done third gemma run"

/general/home/edegreef/gemma-0.98.1-linux-static -bfile $file -k output/relatedness_matrix.cXX.txt -miss 1 -maf 0 -notsnp -bslmm 1 -w 5000000 -s 20000000 -rpace 100 -wpace 1000 -hmin 0 -hmax 1 -rmin 0 -rmax 1 -o $file.bslmm_4
echo "done fourth gemma run"

echo "BSLMMs completed"
