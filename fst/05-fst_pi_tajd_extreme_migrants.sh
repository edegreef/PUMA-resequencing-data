module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0

#fst across genome using extreme 10 early and 10 late migrants, including LD snps and only autosomes. Using non-overlapping 5kb windows
vcftools --vcf /scratch/user/edegreef/reseq/filter_scaffolds/old_vcf_pre-LD/imputed.qualityfiltered.ZWout.recode.vcf --weir-fst-pop early10 --weir-fst-pop late10 --out extreme10.auto.5kb.5kbstep --fst-window-size 5000 --fst-window-step 5000

#preparing vcf for extremes for tajd and pi
vcftools --vcf /scratch/user/edegreef/reseq/filter_scaffolds/old_vcf_pre-LD/imputed.qualityfiltered.ZWout.recode.vcf --keep early10 --recode --recode-INFO-all --out early10.auto
vcftools --vcf /scratch/user/edegreef/reseq/filter_scaffolds/old_vcf_pre-LD/imputed.qualityfiltered.ZWout.recode.vcf --keep late10 --recode --recode-INFO-all --out late10.auto

#tajima's D
vcftools --vcf early10.auto.recode.vcf --TajimaD 5000 --out early10.auto.5kb
vcftools --vcf late10.auto.recode.vcf --TajimaD 5000 --out late10.auto.5kb

#nucleotide diversity
vcftools --vcf early10.auto.recode.vcf --window-pi 5000 --out early10.auto.5kb
vcftools --vcf late10.auto.recode.vcf --window-pi 5000 --out late10.auto.5kb
