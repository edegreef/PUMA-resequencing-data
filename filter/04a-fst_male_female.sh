#!/bin/bash

#Running FST between samples of males and females to identify scaffolds of high differences

input_vcf=imputed.qualityfiltered.recode.vcf
males=list_males
females=list_females
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0

vcftools --vcf $input_vcf --weir-fst-pop $males --weir-fst-pop $females --out males_v_females
