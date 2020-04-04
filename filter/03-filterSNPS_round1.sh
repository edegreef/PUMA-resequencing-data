#!/bin/bash

#Filtering out low-quality and weird SNPs, using kdelmore's scripts as model for most steps
#Inserted a header line using nano command to add "##fileformat=VCFv4.0" to make vcf file readable in more analyses later

#input files:
input_vcf=imputed
r_script=filter.hets.R #for step 5
weird_scaf=weird_scaffolds_RD4_chrom_pos #for step 7

#Step1: Filtering qual > 20 and mapping quality > 20. This filtering step was already done before receiving my "imputed" file but here is some code.
#Need to make sure this line is in header "##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping Quality">" 
#module load vcflib/1.0.0-rc1-GCCcore-6.3.0
#vcffilter -f "QUAL > 20 & MQ > 20" $input_vcf.vcf > $input_vcf.qual.vcf

#Step2: Filtering out indels
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
vcftools --vcf $input_vcf.vcf --remove-indels --recode --recode-INFO-all --out $input_vcf.indels

#Step3: Filtering out minimum allele frequency (maf) < 5%
vcftools --vcf $input_vcf.indels.recode.vcf --maf 0.05 --recode --recode-INFO-all --out $input_vcf.indels.maf

#Step4: Filtering out sites missing genotypes (range 0=missing all, 1=no missing)
vcftools --vcf $input_vcf.indels.maf.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out $input_vcf.indels.maf.0.8miss

#Step5: Filtering out HWE - loci that don't match Hardy-W
vcftools --vcf $input_vcf.indels.maf.0.8miss.recode.vcf --hardy
module load R/3.6.0-iomkl-2018b-recommended-mt
R --vanilla < $r_script
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
vcftools --vcf $input_vcf.indels.maf.0.8miss.recode.vcf --exclude-positions /tmp/tmp_HW_filter --recode --recode-INFO-all --out $input_vcf.indels.maf.0.8miss.het0.6

#Step6: Filtering out triallelic things
awk '{print$1,$2,$5}' $input_vcf.indels.maf.0.8miss.het0.6.recode.vcf > /tmp/tmp1_filter
sed -n '/\,/p' /tmp/tmp1_filter > /tmp/tmp2_filter
sed '/##/d' /tmp/tmp2_filter > /tmp/tmp3_filter
awk -v OFS="\t" '{print$1,$2}' /tmp/tmp3_filter > /tmp/tmp4_filter 
vcftools --vcf $input_vcf.indels.maf.0.8miss.het0.6.recode.vcf --exclude-positions /tmp/tmp4_filter --recode --recode-INFO-all --out $input_vcf.indels.maf.0.8miss.het0.6.tri

#Step7: Filtering out "weird" scaffolds with abnormally high read depths - Using a list determined through other means
vcftools --vcf $input_vcf.indels.maf.0.8miss.het0.6.tri.recode.vcf --exclude-positions weird_scaffolds_RD4_chrom_pos --recode --recode-INFO-all --out $input_vcf.indels.maf.0.8miss.het0.6.tri.weird

#copying/renaming file so it's not so long
cp $input_vcf.indels.maf.0.8miss.het0.6.tri.weird.recode.vcf $input_vcf.qualityfiltered.recode.vcf