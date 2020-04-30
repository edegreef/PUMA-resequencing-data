#!/bin/bash

#Preparing SNPs data and running on smartPCA
snps=/scratch/user/edegreef/reseq/filter_scaffolds/imputed.qualityfiltered.ZWout.pop67.recode.vcf
population=list_4pop_colony_67
samples=list_4pop_indiv_67
output="67_auto"
version="run1_67_auto"

#Converting VCF file to PLINK files (bed, bim, fam)
/general/home/edegreef/plink/plink --vcf $snps --make-bed --allow-extra-chr --set-missing-var-ids @:#\$1,\$2 --out $output

#Editing bim file, changing "scaffold#" to "1", otherwise "bad chrom" error. Can't assign chrom yet because that would remove unlocalized ones.
awk '{print "1",$2,$3,$4,$5,$6}' $output.bim > $output.renamed.bim
sort -k 4 -n $output.renamed.bim > $output.renamed.sort.bim

#Editing fam file to add pop info and duplicating it to last column
awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' $population $output.fam | awk '{print $2, $1, $3, $4, $5, $1}' > $output.popinfo.fam

echo genotypename: $output.bed > genotypename
echo snpname: $output.renamed.sort.bim > snpname
echo indivname: $output.popinfo.fam > indivname
echo evecoutname: smartpca.$version.evec.txt > evecoutname
echo evaloutname: smartpca.$version.eval.txt > evaloutname

cat genotypename snpname indivname evecoutname evaloutname > pcaparfile_$version
rm genotypename snpname indivname evecoutname evaloutname

#use the pcarparfile for running smartpca
module load EIGENSOFT/7.2.1-intel-2019a
smartpca -p pcaparfile_$version > smartpca.log

mv smartpca.log smartpca.$version.log

#use the evec file for plotting in R
