#library(BiocManager)
#BiocManager::install("LEA")

library(LEA)

#convert vcf to lfmm format
vcf2lfmm("Imputed_all_snps4615223.77.recode.vcf") #automatically puts output in working directory

project=NULL
project=lfmm("Imputed_all_snps4615223.77.recode.lfmm", "bloom77.cov5.col.year.pc1.sex.age.env", K=1, repetitions=5, project="new", all=TRUE)
p=lfmm.pvalues(project, K=1)

save(p, file="bloom77.cov5.RData")