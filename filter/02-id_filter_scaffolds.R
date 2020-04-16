library(tidyverse)

######## Setting up input files
setwd("E:/Evelien's Dropbox/Dropbox/PUMA/Bioinformatics/scaffold_filter")

#list of all scaffolds & positions of vcf file
list_all <- read.table("list_all", header=T)

#Scaffold read depths between males and females were obtained using raw fastq files and reference genome
data <- read.csv("scaffold_read_depths.csv")
data <- select(data,
               scaffold,
               length,
               avg_female_RD,
               avg_male_RD,
               male_female=m.f.m.f, #this is male-female/male+female
               avg_RD)



######## Identifying "weird" scaffolds with high read depths
#selecting scaffolds based on average read depths
weird_scaffolds <- filter(data, avg_RD >= 4)
sum(weird_scaffolds$length) #total length in bp

#extract scaffold column
weird_scaf <- select(weird_scaffolds, scaffold)

#filtering list_all based on scaffolds and merging columns
filtered <- subset(list_all, CHROM %in% weird_scaf$scaffold)
filtered2 <- paste(filtered$CHROM, filtered$POS)



######## Identifying Z-linked scaffolds
#males have 2x read depths than females for Z chromosome. Using parameter in "male_female" column representing male-female/male+female.
weird_removed <- filter(data, avg_RD < 4)
Z_scaffolds <- filter(weird_removed, male_female > 0.04)
sum(Z_scaffolds$length) #total length in bp

#extract scaffold column
Z_scaf <- select(Z_scaffolds, scaffold)

#filtering list_all based on scaffolds and merging columns
filtered_Z <- subset(list_all, CHROM %in% Z_scaf$scaffold)
filtered_Z2 <- paste(filtered_Z$CHROM, filtered_Z$POS)



######## Identifying W-linked scaffolds
#females have a W chromosome and males do not. Using paramters of male and female read depths to filter scaffolds out.
W_scaffolds <- filter(weird_removed, avg_male_RD < 0.7 & avg_female_RD > 0.8)
sum(W_scaffolds$length) #total length in bp

#extract scaffold column
W_scaf <- select(W_scaffolds, scaffold)

#filtering list_all based on scaffolds and merging columns
filtered_W <- subset(list_all, CHROM %in% W_scaf$scaffold)
filtered_W2 <- paste(filtered_W$CHROM, filtered_W$POS)



######## Saving list files
#using these lists later in vcftools
write.csv(filtered2,"weird_scaffolds_RD4_chrom_pos.csv", row.names=FALSE)
write.csv(filtered_Z2,"Z_scaffolds_MF0.04_chrom_pos.csv", row.names=FALSE)
write.csv(filtered_W2,"W_scaffolds_M0.7_F0.8_chrom_pos.csv", row.names=FALSE)

#saving these for records (includes read depths)
write.csv(weird_scaffolds,'weird_scaffolds_RD4.csv')
write.csv(Z_scaffolds,'Z_scaffolds_MF0.04.csv')
write.csv(W_scaffolds,'W_scaffolds_M0.7_F0.8.csv')
