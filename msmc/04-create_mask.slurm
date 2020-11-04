#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=mask_file
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

ref_genome=/scratch/user/edegreef/ref_genome/PUMA2_03022020.fasta
input_bam=/scratch/user/edegreef/reseq_bam/extra_subsp_bams/PM150_S462.duprem.bam
ID=PM150_min1MB
list=list_scaf_min1MB
coverage=14

module load SAMtools/1.10-iccifort-2019.5.281
module load BCFtools/1.10.2-iccifort-2019.5.281
module load Python/3.8.2-GCCcore-9.3.0

cd /scratch/user/edegreef/msmc

mkdir mask_output_$ID

while read scaffold
do
	samtools mpileup -q 20 -Q 20 -C 50 -u -r $scaffold -f $ref_genome $input_bam | bcftools call -c -V indels | msmc-tools/bamCaller.py $coverage mask_output_$ID/$ID.mask.$scaffold.bed.gz | gzip -c > mask_output_$ID/$ID.$scaffold.vcf.gz
done < $list




