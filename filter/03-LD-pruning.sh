#!/bin/bash

#filtering out SNPs that are in linkage disequilibrium (LD pruning)

#ID snps in LD
/general/home/edegreef/plink/plink --allow-extra-chr --indep-pairwise 200 20 0.2 --out imputed.qualityfiltered --set-missing-var-ids @:#\$1,\$2 --vcf imputed.qualityfiltered.recode.vcf

#editing prune.in list to run in vcftools
sed 's/:/\t/g' imputed.qualityfiltered.prune.in > temp.in
sed 's/...$//' temp.in > temp2.in
rm imputed.qualityfiltered.prune.in
mv temp2.in imputed.qualityfiltered.prune.in

#LD pruning snps!
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
vcftools --vcf imputed.qualityfiltered.recode.vcf --positions imputed.qualityfiltered.prune.in --recode --recode-INFO-all --out imputed.qualityfiltered.LD
