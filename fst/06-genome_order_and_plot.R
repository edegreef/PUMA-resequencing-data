#plotting fst, tajd, pi across genome

setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/Bioinformatics/reseq/other_stats")

#load genome chromosome order info
scafs <- read.csv("purple_scafs_ordered_by_chicken.csv") #scaffold/chromosome info
scafs$id <- 1:nrow(scafs)
scafs$purple_scaf <- gsub("Q","",scafs$purple_scaf)

#load fst data file
fst <- read.table("FST/ABFL.auto.5kb.5kbstep.windowed.weir.fst", header=TRUE)

fst$best_chick_chr <- as.character(NA)
fst$best_chick_id <- as.numeric(0)
fst$orientation <- as.character(NA)

#adding chr by matching scaffold name
for (i in 1:length(fst$CHROM)){
  for (j in 1:length(scafs$chr_num)){
    
    if (fst$CHROM[i] == scafs$purple_scaf[j]){
      fst$best_chick_chr[i] <- as.character(scafs$chr_num[j])
      fst$best_chick_id[i] <- as.numeric(scafs$id[j])
      fst$orientation[i] <- as.character(scafs$ori[j])
    }
  }
}

write.csv(fst,'FST/ABFL.auto.5kb.5kbstep.windowed.weir.fst.orderedchr.csv') #save output as backup

summary(fst$best_chick_chr)
fst$new_chr <- as.numeric(as.character(factor(fst$best_chick_chr,
                                              levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '33', 'W', 'Z'), labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)))) 
#W=30
#Z=31
#U=32

fst$pos <- ((fst$BIN_START + fst$BIN_END) /2)
fst$new_chr[is.na(fst$new_chr)] <- 32

#sorting by new_chr, then by id
fst_order <- fst[order(fst[,10], fst[,8]),]
data_ordered <- fst_order

#need to order by chrom & position in each chr
chr1 <- subset(data_ordered, new_chr==1)
chr1$order_number <- 1:nrow(chr1)
chr2 <- subset(data_ordered, new_chr==2)
chr2$order_number <- 1:nrow(chr2)
chr3 <- subset(data_ordered, new_chr==3)
chr3$order_number <- 1:nrow(chr3)
chr4 <- subset(data_ordered, new_chr==4)
chr4$order_number <- 1:nrow(chr4)
chr5 <- subset(data_ordered, new_chr==5)
chr5$order_number <- 1:nrow(chr5)
chr6 <- subset(data_ordered, new_chr==6)
chr6$order_number <- 1:nrow(chr6)
chr7 <- subset(data_ordered, new_chr==7)
chr7$order_number <- 1:nrow(chr7)
chr8 <- subset(data_ordered, new_chr==8)
chr8$order_number <- 1:nrow(chr8)
chr9 <- subset(data_ordered, new_chr==9)
chr9$order_number <- 1:nrow(chr9)
chr10 <- subset(data_ordered, new_chr==10)
chr10$order_number <- 1:nrow(chr10)
chr11 <- subset(data_ordered, new_chr==11)
chr11$order_number <- 1:nrow(chr11)
chr12 <- subset(data_ordered, new_chr==12)
chr12$order_number <- 1:nrow(chr12)
chr13 <- subset(data_ordered, new_chr==13)
chr13$order_number <- 1:nrow(chr13)
chr14 <- subset(data_ordered, new_chr==14)
chr14$order_number <- 1:nrow(chr14)
chr15 <- subset(data_ordered, new_chr==15)
chr15$order_number <- 1:nrow(chr15)
chr16 <- subset(data_ordered, new_chr==16)
chr16$order_number <- 1:nrow(chr16)
chr17 <- subset(data_ordered, new_chr==17)
chr17$order_number <- 1:nrow(chr17)
chr18 <- subset(data_ordered, new_chr==18)
chr18$order_number <- 1:nrow(chr18)
chr19 <- subset(data_ordered, new_chr==19)
chr19$order_number <- 1:nrow(chr19)
chr20 <- subset(data_ordered, new_chr==20)
chr20$order_number <- 1:nrow(chr20)
chr21 <- subset(data_ordered, new_chr==21)
chr21$order_number <- 1:nrow(chr21)
chr22 <- subset(data_ordered, new_chr==22)
chr22$order_number <- 1:nrow(chr22)
chr23 <- subset(data_ordered, new_chr==23)
chr23$order_number <- 1:nrow(chr23)
chr24 <- subset(data_ordered, new_chr==24)
chr24$order_number <- 1:nrow(chr24)
chr25 <- subset(data_ordered, new_chr==25)
chr25$order_number <- 1:nrow(chr25)
chr26 <- subset(data_ordered, new_chr==26)
chr26$order_number <- 1:nrow(chr26)
chr27 <- subset(data_ordered, new_chr==27)
chr27$order_number <- 1:nrow(chr27)
chr28 <- subset(data_ordered, new_chr==28)
chr28$order_number <- 1:nrow(chr28)
chr29 <- subset(data_ordered, new_chr==29)
chr29$order_number <- 1:nrow(chr29)
chr30 <- subset(data_ordered, new_chr==30)
chr30$order_number <- 1:nrow(chr30)
chr31 <- subset(data_ordered, new_chr==31)
chr31$order_number <- 1:nrow(chr31)
chr32 <- subset(data_ordered, new_chr==32)
chr32$order_number <- 1:nrow(chr32)

chr_combined <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chr26, chr27, chr28, chr29)

chr_combined$plot_order <- 1:nrow(chr_combined)

write.csv(chr_combined, "FST/ABFL.auto.5kb.5kbstep.windowed.weir.fst.orderedchr.plotorder.GOOD.csv")
for_plot <- chr_combined

#check to make sure chr and pos order is good
library(qqman)
png("ABFfst window 5kb.png", w=1200, h=400)
manhattan(for_plot,main="FST extreme10 5kb window test", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="CHROM", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), cex.axis=0.9, suggestiveline = F, genomewideline = F)
dev.off()

#next going to subtract AB vs FL fst from early vs late fst
ABFL <- chr_combined
EX10 <- read.csv("C:/Users/Evelien de Greef/Dropbox/PUMA/Bioinformatics/reseq/other_stats/FST/extreme10.auto.5kb.5kbstep.windowed.weir.fst.orderedchr.plotorder.GOOD.csv", header=T)

ABFL$snp <- paste(ABFL$CHROM, ABFL$BIN_START)
EX10$snp <- paste(EX10$CHROM, EX10$BIN_START)

library(dplyr)
ABFL_subset <- ABFL %>%
  filter(snp %in% EX10$snp)

EX10_subset <- EX10 %>%
  filter(snp %in% ABFL$snp)

EX10_subset$ABFL_fst <- ABFL_subset$WEIGHTED_FST
EX10_subset$fst_diff <- EX10_subset$WEIGHTED_FST - EX10_subset$ABFL_fst

png("C:/Users/Evelien de Greef/Dropbox/PUMA/Bioinformatics/reseq/other_stats/FST/extreme10_ABFL_fstdifference.png", w=1200, h=400)
manhattan(EX10_subset,main="difference between extreme10 and ABFL fst", chr="new_chr", bp="order_number", p="fst_diff", snp="CHROM", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), cex.axis=0.9, suggestiveline = F, genomewideline = F)
dev.off()


#some code for zooming in at particular spots

highlightsnps <- as.character(c(175650:175730)) #+/- 40bp

#zooming in on chr 24
manhattan(subset(for_plot, new_chr==24),main="extreme10 FST S129:449783", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7),cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight=highlightsnps)

chr1 <- subset(for_plot, new_chr==1)
chr4 <- subset(for_plot, new_chr==4)
chr2 <- subset(for_plot, new_chr==2)

orderfst<- chr1[order(chr1$WEIGHTED_FST),]
#looking at the fsts > 0.3

plot <- subset(chr1, plot_order > 5500 & plot_order < 6200)
plot(plot$plot_order, plot$WEIGHTED_FST)
abline(v=c(5714, 5718), col="blue")
abline(v=c(5731), col="orange")
abline(v=c(5791, 5792), col="blue")
abline(v=c(5863, 5865), col="orange")
abline(v=c(5866,5867), col="blue")
abline(v=c(5903, 5915), col="orange")
abline(v=c(5918, 5919), col="blue")
abline(v=c(5920, 5928), col="orange")
abline(v=c(5938, 5948), col="blue")
abline(v=c(5953, 5954), col="orange")
abline(v=c(5977), col="blue")
abline(v=c(5985, 5988), col="orange")
abline(v=c(5992, 5994), col="blue")
abline(v=c(6141, 6145), col="orange")
abline(v=c(6146, 6147), col="blue")

##################################################
####### trying out zooming in on other snps
whole_fst <- read.csv("extreme20.auto.5kb.5kbstep.windowed.weir.fst.orderedchr.plotorder.GOOD.csv", header=TRUE)

manhattan(whole_fst,main="test", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="X", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), cex.axis=0.9, suggestiveline = F, genomewideline = F)

##### snp 1
chr24to26 <- subset(whole_fst, new_chr >= 24 & new_chr <= 26)
highlightsnps <- as.character(c(176551:176630)) #20bp left and right

png("S214-738983 - extreme20 5kb.png", w=1200, h=400)
manhattan(chr24to26, main="S214:738983 - extreme20 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)
dev.off()

##### snp 2
scaf13 <- subset(whole_fst, CHROM=="scaffold13")
highlightsnps <- as.character(c(77484:77564)) #40bp left and right

png("S13-2452263 - extreme20 5kb.png", w=1200, h=400)
manhattan(subset(whole_fst, new_chr==3), main="S13:2452263 - extreme20 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), xlim=c(10000, 15000), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)
dev.off()

##### snp 3
highlightsnps <- as.character(c(110471:110551)) #40bp left and right

png("S39-10682289 - extreme20 5kb.png", w=1200, h=400)
manhattan(subset(whole_fst, new_chr==5), main="S39:10682289 - extreme20 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), xlim=c(2000, 10000), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)
dev.off()


##### snp 4
highlightsnps <- as.character(c(140678:140748)) #40bp left and right

png("S31-3217866 - extreme20 5kb.png", w=1200, h=400)
manhattan(subset(whole_fst, new_chr==9), main="S31:3217866 - extreme20 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7),xlim=c(500, 3500), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)
dev.off()

##### snp 5
highlightsnps <- as.character(c(104886:104966)) #40bp left and right

png("S9-3837598 - extreme20 5kb.png", w=1200, h=400)
manhattan(subset(whole_fst, new_chr==5), main="S9:3837598 - extreme20 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7),xlim=c(0,4000), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)
dev.off()

##### snp 6
highlightsnps <- as.character(c(30530:30610)) #40bp left and right
png("S2-18994041 - extreme20 5kb.png", w=1200, h=400)
manhattan(subset(whole_fst, new_chr==1), main="S2:18994041 - extreme20 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), xlim=c(28000,32000), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)
dev.off()


##### snp 7
highlightsnps <- as.character(c(43228:43308)) #40bp left and right
png("S10-31178211 - extreme20 5kb.png", w=1200, h=400)
manhattan(subset(whole_fst, new_chr==2), main="S10:31178211 - extreme20 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), xlim=c(2500,10000), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)
dev.off()


##### snp 8
highlightsnps <- as.character(c(102876:102956)) #40bp left and right
png("S3-15038051 - extreme20 5kb.png", w=1200, h=400)
manhattan(subset(whole_fst, new_chr==4), main="S3:15038051 - extreme20 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), xlim=c(13000, 16000),cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)
dev.off()

#### snp12:405479
chr4 <- subset(for_plot, new_chr==4)
highlightsnps <- as.character(c(93290:93332))
manhattan(chr4, main="S12:405479", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.5),xlim=c(4500,6500),cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)


#### snp37:4184869
chr2 <- subset(for_plot, new_chr==2)
highlightsnps <- as.character(c(53515:53555))
manhattan(chr2, main="S37:4184869", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.5),xlim=c(17000,18200),cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)



#going to try highlighting all snps on the plot
highlightsnps <- as.character(c(176590, 77524, 110511, 140708, 104926, 30570, 43268, 102916)) 
#40 positions on each side
highlightsnps_chunk <- as.character(c(176550:176630, 
                                      77484:77564, 
                                      110471:110551, 
                                      140668:140748,
                                      104886:104966,
                                      30530:30610,
                                      43228:43308,
                                      102876:102956))

png("allsnp- extreme10 5kb.png", w=1200, h=400)
manhattan(whole_fst,main="all - extreme10 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps)
dev.off()

png("allsnpchunks- extreme20 5kb.png", w=1200, h=400)
manhattan(whole_fst,main="all - extreme20 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps_chunk)
dev.off()


png("allsnpchunks- extreme10 5kb - large.png", w=4000, h=1200, point=50)
manhattan(whole_fst,main="all - extreme10 5kb", chr="new_chr", bp="order_number", p="WEIGHTED_FST", snp="plot_order", logp=FALSE, ylab="Weir and Cockerham",  ylim = c(-0.05, 0.7), cex.axis=0.9, suggestiveline = F, genomewideline = F, highlight = highlightsnps_chunk)
dev.off()

########################################################################
########################################################################
#tajima's D

scafs <- read.csv("purple_scafs_ordered_by_chicken.csv") #scaffold/chromosome info
scafs$id <- 1:nrow(scafs)
scafs$purple_scaf <- gsub("Q","",scafs$purple_scaf)

TD <- read.table('early10.auto.5kb.Tajima.D', header=T)

TD$best_chick_chr <- as.character(NA)
TD$best_chick_id <- as.numeric(0)
TD$orientation <- as.character(NA)

for (i in 1:length(TD$CHROM)){
  for (j in 1:length(scafs$chr_num)){
    
    if (TD$CHROM[i] == scafs$purple_scaf[j]){
      TD$best_chick_chr[i] <- as.character(scafs$chr_num[j])
      TD$best_chick_id[i] <- as.numeric(scafs$id[j])
      TD$orientation[i] <- as.character(scafs$ori[j])
    }
  }
}

write.csv(TD,'early10.auto.5kb.Tajima.D.csv')

###### loading ordered csv
TD <- read.csv('late10.auto.5kb.Tajima.D.orderedchr.csv', header=T)

summary(TD$best_chick_chr)
TD$new_chr <- as.numeric(as.character(factor(TD$best_chick_chr,
                                              levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '33', 'W', 'Z'), labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)))) 
#W=30
#Z=31
#U=32

TD$pos <- ((TD$BIN_START + 10000) /2)
TD$new_chr[is.na(TD$new_chr)] <- 32

#sorting by new_chr, then by id
TD <- subset( TD, select = -X )

TD_order <- TD[order(TD[,8], TD[,6]),]
TD_order$order_number <- 1:nrow(TD_order)
data_ordered <- TD_order

#need to order by chrom & position in each chr
chr1 <- subset(data_ordered, new_chr==1)
chr1$order_number <- 1:nrow(chr1)
chr2 <- subset(data_ordered, new_chr==2)
chr2$order_number <- 1:nrow(chr2)
chr3 <- subset(data_ordered, new_chr==3)
chr3$order_number <- 1:nrow(chr3)
chr4 <- subset(data_ordered, new_chr==4)
chr4$order_number <- 1:nrow(chr4)
chr5 <- subset(data_ordered, new_chr==5)
chr5$order_number <- 1:nrow(chr5)
chr6 <- subset(data_ordered, new_chr==6)
chr6$order_number <- 1:nrow(chr6)
chr7 <- subset(data_ordered, new_chr==7)
chr7$order_number <- 1:nrow(chr7)
chr8 <- subset(data_ordered, new_chr==8)
chr8$order_number <- 1:nrow(chr8)
chr9 <- subset(data_ordered, new_chr==9)
chr9$order_number <- 1:nrow(chr9)
chr10 <- subset(data_ordered, new_chr==10)
chr10$order_number <- 1:nrow(chr10)
chr11 <- subset(data_ordered, new_chr==11)
chr11$order_number <- 1:nrow(chr11)
chr12 <- subset(data_ordered, new_chr==12)
chr12$order_number <- 1:nrow(chr12)
chr13 <- subset(data_ordered, new_chr==13)
chr13$order_number <- 1:nrow(chr13)
chr14 <- subset(data_ordered, new_chr==14)
chr14$order_number <- 1:nrow(chr14)
chr15 <- subset(data_ordered, new_chr==15)
chr15$order_number <- 1:nrow(chr15)
chr16 <- subset(data_ordered, new_chr==16)
chr16$order_number <- 1:nrow(chr16)
chr17 <- subset(data_ordered, new_chr==17)
chr17$order_number <- 1:nrow(chr17)
chr18 <- subset(data_ordered, new_chr==18)
chr18$order_number <- 1:nrow(chr18)
chr19 <- subset(data_ordered, new_chr==19)
chr19$order_number <- 1:nrow(chr19)
chr20 <- subset(data_ordered, new_chr==20)
chr20$order_number <- 1:nrow(chr20)
chr21 <- subset(data_ordered, new_chr==21)
chr21$order_number <- 1:nrow(chr21)
chr22 <- subset(data_ordered, new_chr==22)
chr22$order_number <- 1:nrow(chr22)
chr23 <- subset(data_ordered, new_chr==23)
chr23$order_number <- 1:nrow(chr23)
chr24 <- subset(data_ordered, new_chr==24)
chr24$order_number <- 1:nrow(chr24)
chr25 <- subset(data_ordered, new_chr==25)
chr25$order_number <- 1:nrow(chr25)
chr26 <- subset(data_ordered, new_chr==26)
chr26$order_number <- 1:nrow(chr26)
chr27 <- subset(data_ordered, new_chr==27)
chr27$order_number <- 1:nrow(chr27)
chr28 <- subset(data_ordered, new_chr==28)
chr28$order_number <- 1:nrow(chr28)
chr29 <- subset(data_ordered, new_chr==29)
chr29$order_number <- 1:nrow(chr29)
chr30 <- subset(data_ordered, new_chr==30)
chr30$order_number <- 1:nrow(chr30)
chr31 <- subset(data_ordered, new_chr==31)
chr31$order_number <- 1:nrow(chr31)
chr32 <- subset(data_ordered, new_chr==32)
chr32$order_number <- 1:nrow(chr32)

chr_combined <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chr26, chr27, chr28, chr29)

chr_combined$plot_order <- 1:nrow(chr_combined)

write.csv(chr_combined, "late10.auto.5kb.Tajima.D.orderedchr.plotorder.good.csv")

library(ggplot2)

ggplot(chr_combined, aes(x=plot_order, y=TajimaD))+
  geom_point(color="gray", alpha=0.7) +
  theme_classic()

manhattan(chr_combined,main="TD-late - extreme10 5kb", chr="new_chr", bp="order_number", p="TajimaD", snp="plot_order", logp=FALSE, ylab="TajimaD",  ylim = c(-3,3), cex.axis=0.9, suggestiveline = F, genomewideline = F)

setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/Bioinformatics/reseq/other_stats/TAJD")
TD_early <- read.csv("early10.auto.5kb.Tajima.D.orderedchr.plotorder.good.csv", header=T)
TD_late <- read.csv("late10.auto.5kb.Tajima.D.orderedchr.plotorder.good.csv", header=T)



########################################################################
########################################################################
#now same thing with pi

pi <- read.table('early10.auto.5kb.windowed.pi', header=T)

pi$best_chick_chr <- as.character(NA)
pi$best_chick_id <- as.numeric(0)
pi$orientation <- as.character(NA)

for (i in 1:length(pi$CHROM)){
  for (j in 1:length(scafs$chr_num)){
    
    if (pi$CHROM[i] == scafs$purple_scaf[j]){
      pi$best_chick_chr[i] <- as.character(scafs$chr_num[j])
      pi$best_chick_id[i] <- as.numeric(scafs$id[j])
      pi$orientation[i] <- as.character(scafs$ori[j])
    }
  }
}


write.csv(pi,'late10.auto.5kb.windowed.pi.orderedchr.csv')

pi <- read.csv('late10.auto.5kb.windowed.pi.orderedchr.csv', header=T)

summary(pi$best_chick_chr)
pi$new_chr <- as.numeric(as.character(factor(pi$best_chick_chr,
                                             levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '33', 'W', 'Z'), labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)))) 
#W=30
#Z=31
#U=32

pi <- subset( pi, select = -X )

pi$pos <- ((pi$BIN_START + pi$BIN_END) /2)
pi$new_chr[is.na(pi$new_chr)] <- 32

#sorting by new_chr, then by id
pi_order <- pi[order(pi[,9], pi[,7]),]
pi_order$order_number <- 1:nrow(pi_order)
data_ordered <- pi_order

#need to order by chrom & position in each chr
chr1 <- subset(data_ordered, new_chr==1)
chr1$order_number <- 1:nrow(chr1)
chr2 <- subset(data_ordered, new_chr==2)
chr2$order_number <- 1:nrow(chr2)
chr3 <- subset(data_ordered, new_chr==3)
chr3$order_number <- 1:nrow(chr3)
chr4 <- subset(data_ordered, new_chr==4)
chr4$order_number <- 1:nrow(chr4)
chr5 <- subset(data_ordered, new_chr==5)
chr5$order_number <- 1:nrow(chr5)
chr6 <- subset(data_ordered, new_chr==6)
chr6$order_number <- 1:nrow(chr6)
chr7 <- subset(data_ordered, new_chr==7)
chr7$order_number <- 1:nrow(chr7)
chr8 <- subset(data_ordered, new_chr==8)
chr8$order_number <- 1:nrow(chr8)
chr9 <- subset(data_ordered, new_chr==9)
chr9$order_number <- 1:nrow(chr9)
chr10 <- subset(data_ordered, new_chr==10)
chr10$order_number <- 1:nrow(chr10)
chr11 <- subset(data_ordered, new_chr==11)
chr11$order_number <- 1:nrow(chr11)
chr12 <- subset(data_ordered, new_chr==12)
chr12$order_number <- 1:nrow(chr12)
chr13 <- subset(data_ordered, new_chr==13)
chr13$order_number <- 1:nrow(chr13)
chr14 <- subset(data_ordered, new_chr==14)
chr14$order_number <- 1:nrow(chr14)
chr15 <- subset(data_ordered, new_chr==15)
chr15$order_number <- 1:nrow(chr15)
chr16 <- subset(data_ordered, new_chr==16)
chr16$order_number <- 1:nrow(chr16)
chr17 <- subset(data_ordered, new_chr==17)
chr17$order_number <- 1:nrow(chr17)
chr18 <- subset(data_ordered, new_chr==18)
chr18$order_number <- 1:nrow(chr18)
chr19 <- subset(data_ordered, new_chr==19)
chr19$order_number <- 1:nrow(chr19)
chr20 <- subset(data_ordered, new_chr==20)
chr20$order_number <- 1:nrow(chr20)
chr21 <- subset(data_ordered, new_chr==21)
chr21$order_number <- 1:nrow(chr21)
chr22 <- subset(data_ordered, new_chr==22)
chr22$order_number <- 1:nrow(chr22)
chr23 <- subset(data_ordered, new_chr==23)
chr23$order_number <- 1:nrow(chr23)
chr24 <- subset(data_ordered, new_chr==24)
chr24$order_number <- 1:nrow(chr24)
chr25 <- subset(data_ordered, new_chr==25)
chr25$order_number <- 1:nrow(chr25)
chr26 <- subset(data_ordered, new_chr==26)
chr26$order_number <- 1:nrow(chr26)
chr27 <- subset(data_ordered, new_chr==27)
chr27$order_number <- 1:nrow(chr27)
chr28 <- subset(data_ordered, new_chr==28)
chr28$order_number <- 1:nrow(chr28)
chr29 <- subset(data_ordered, new_chr==29)
chr29$order_number <- 1:nrow(chr29)
chr30 <- subset(data_ordered, new_chr==30)
chr30$order_number <- 1:nrow(chr30)
chr31 <- subset(data_ordered, new_chr==31)
chr31$order_number <- 1:nrow(chr31)
chr32 <- subset(data_ordered, new_chr==32)
chr32$order_number <- 1:nrow(chr32)

chr_combined <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chr26, chr27, chr28, chr29)

chr_combined$plot_order <- 1:nrow(chr_combined)

write.csv(chr_combined, "late10.auto.5kb.windowed.pi.orderedchr.plotorder.good.csv")

library(ggplot2)

chr_combined <- read.csv("early10.auto.5kb.windowed.pi.orderedchr.plotorder.good.csv", header=T)

ggplot(chr_combined, aes(x=plot_order, y=PI))+
  geom_point(color="gray", alpha=0.5) +
  theme_classic()

png("PI_early10_5kb_1500x500.png", width=1500, height=500)
manhattan(chr_combined,main="pi early - extreme10 5kb", chr="new_chr", bp="order_number", p="PI", snp="plot_order", logp=FALSE, ylab="PI",  ylim = c(-0.001,0.015), cex.axis=0.9, suggestiveline = F, genomewideline = F)
dev.off()


#############################################
#trying density addition, just playing around with stuff
library(MASS)
library(viridis)

theme_set(theme_bw(base_size=16))
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

set.seed(1)
dat <- data.frame(
  x = c(
    rnorm(1e4, mean = 0, sd = 0.1),
    rnorm(1e3, mean = 0, sd = 0.1)
  ),
  y = c(
    rnorm(1e4, mean = 0, sd = 0.1),
    rnorm(1e3, mean = 0.1, sd = 0.2)
  )
)

ggplot(chr_combined) + geom_point(aes(plot_order, PI))

chr_combined$density <- get_density(chr_combined$plot_order, chr_combined$PI, n = 100)
ggplot(chr_combined) + geom_point(aes(plot_order, PI, color = density)) + scale_color_viridis()

#trying highlighting specific points
#first going to subset data with chr <=29
chr_part <- subset(chr_combined, new_chr <=29)
ggplot(chr_part, aes(x=plot_order, y=PI))+
  geom_point(color="gray", alpha=0.5) +
  theme_classic()



chr_part %>%
  ggplot(aes(x=plot_order, y=PI))+
  geom_point()

library(tidyverse)
chr_part %>%
  mutate(highlight=ifelse(N_VARIANTS > 300, T,F)) %>%
  ggplot(aes(x=plot_order, y=PI))+
  geom_point(aes(color=highlight))


TDtest <- read.csv("window_100kb/early20.allchrom.100kb.LD.Tajima.D.orderedchr.plotorder.csv", header=T)
TD_part <- subset(TDtest, new_chr <=29)
TD_part %>%
  mutate(highlight=ifelse(N_SNPS > 300, T,F)) %>%
  ggplot(aes(x=plot_order, y=TajimaD))+
  geom_point(aes(color=highlight))
