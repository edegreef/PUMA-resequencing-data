#!/bin/bash

#Script for getting all FST values among the 4 Purple Martin populations (AB, PA, TX, FL)
#Need to include lists per population containing sample IDs.

input_vcf=/scratch/user/edegreef/reseq/filter_scaffolds/imputed.qualityfiltered.ZWout.pop.recode.vcf #autosomes
type=auto

module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0

#AB vs FL
vcftools --vcf $input_vcf --weir-fst-pop alberta --weir-fst-pop florida --out AB_vs_FL_$type --fst-window-size 5000 --fst-window-step 5000

#AB vs PA
vcftools --vcf $input_vcf --weir-fst-pop alberta --weir-fst-pop pennsylvania --out AB_vs_PA_$type --fst-window-size 5000 --fst-window-step 5000

#AB vs TX
vcftools --vcf $input_vcf --weir-fst-pop alberta --weir-fst-pop texas --out AB_vs_TX_$type --fst-window-size 5000 --fst-window-step 5000

#FL vs PA
vcftools --vcf $input_vcf --weir-fst-pop florida --weir-fst-pop pennsylvania --out FL_vs_PA_$type --fst-window-size 5000 --fst-window-step 5000

#FL vs TX
vcftools --vcf $input_vcf --weir-fst-pop florida --weir-fst-pop texas --out FL_vs_TX_$type --fst-window-size 5000 --fst-window-step 5000

#PA vs TX
vcftools --vcf $input_vcf --weir-fst-pop pennsylvania --weir-fst-pop texas --out PA_vs_TX_$type --fst-window-size 5000 --fst-window-step 5000

echo "all FST runs done"