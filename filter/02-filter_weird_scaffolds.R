setwd("E:/Evelien's Dropbox/Dropbox/PUMA/Bioinformatics/scaffold_filter")

#Scaffold read depths from samples were obtained using raw bam files and reference genome
data <- read.csv("scaffold_read_depths.csv")

library(tidyverse)

#renaming column names
data <- select(data,
               scaffold,
               length,
               avg_female_RD,
               avg_male_RD,
               male_female=m.f.m.f, #this is male-female/male+female
               avg_RD)

#making a list of "weird" scaffolds with average RD >= 4
weird_scaffolds <- filter(data, avg_RD >= 4)

#just to see how many bp
sum(weird_scaffolds$length)

#save for record
write.csv(weird_scaffolds,'weird_scaffolds_RD4.csv')

#extract scaffold column
scaffolds <- select(weird_scaffolds, scaffold)

#loading list of all scaffolds & positions of vcf file
list_all <- read.table("list_all", header=T)

#filtering list_all based on scaffolds and merging columns
filtered <- subset(list_all, CHROM %in% scaffolds$scaffold)
filtered2 <- paste(filtered$CHROM, filtered$POS)
write.csv(filtered2,"weird_scaffolds_RD4_chrom_pos.csv", row.names=FALSE)