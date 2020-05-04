#!/usr/bin/Rscript

#script for plotting FST by chromosome order
setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/Bioinformatics/reseq/FST")

fst <- read.table("FL_vs_TX_auto.67.LD.windowed.weir.fst", header=T) #fst data
scafs <- read.csv("/Users/Evelien de Greef/Dropbox/PUMA/Bioinformatics/reseq/gemma/purple_scafs_ordered_by_chicken.csv") #scaffold/chr info
scafs$id <- 1:nrow(scafs)
scafs$purple_scaf <- gsub("Q","",scafs$purple_scaf)

#making new columns for adding chromosome info
fst$best_chick_chr <- as.character(NA)
fst$best_chick_id <- as.numeric(0)
fst$orientation <- as.character(NA)

#adding chr by matching scaffold name, this next chunk takes a while to run
for (i in 1:length(fst$CHROM)){
  for (j in 1:length(scafs$chr_num)){
    
    if (fst$CHROM[i] == scafs$purple_scaf[j]){
      fst$best_chick_chr[i] <- as.character(scafs$chr_num[j])
      fst$best_chick_id[i] <- as.numeric(scafs$id[j])
      fst$orientation[i] <- as.character(scafs$ori[j])
    }
  }
}

write.csv(fst,'FL_vs_TX_auto.67.LD.windowed.weir.fst.orderedchr.csv')




######### after ordering chromosomes
fst$pos <- ((fst$BIN_START + fst$BIN_END) /2) #this step necessary for windowed fsts with bin start and end
fst$best_chick_chr[is.na(fst$best_chick_chr)] <- 34
fst$best_chick_chr <- as.numeric(fst$best_chick_chr)

#quick plot weighted fst
plot(fst$pos, fst$WEIGHTED_FST, main="AB/FL windowed FST")

#manhattan plot
library(qqman)

manhattan(fst,main="AB vs PA windowed FST", chr="best_chick_chr", bp="pos", p="WEIGHTED_FST", snp="CHROM", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.1, 1), cex.axis=0.9, suggestiveline = F, genomewideline = F)
