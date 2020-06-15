setwd("E:/Evelien's Dropbox/Dropbox/PUMA/Bioinformatics/scaffold_filter")

fst <- read.table("males_v_females.weir.fst", header=T)

#just looking at some things
hist(fst$WEIR_AND_COCKERHAM_FST)
fst_sort <- fst[order(fst$WEIR_AND_COCKERHAM_FST),]

#filtering out scaffold#s that may be sex-linked (fst >= 0.2)
library(tidyverse)
sex_linked <- filter(fst, WEIR_AND_COCKERHAM_FST >= 0.2)

sex_linked_chrom_pos <- paste(sex_linked$CHROM, sex_linked$POS)
write.csv(sex_linked_chrom_pos,'sex_linked_chrom_pos.csv')
