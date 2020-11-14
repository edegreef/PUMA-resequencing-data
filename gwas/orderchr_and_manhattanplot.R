#!/usr/bin/Rscript

#ordering bslmms and lmms for manhattan plots

#### transferring over scaffold/chromosome info to dataframe

snplist <- read.csv("imputed.qualityfiltered.snplist.CHROM.POS.orderchr.csv") #load snp list with chrom info
data <- read.csv("geobirds87.ZWoutall.springresid_para.csv", header=T) 
id="geobirds87.ZWoutall.springresid_para"

#make snp column in snplist and data with unique snp ids
data$snp <- paste(data$CHR, data$PS)
snplist$snp <- paste(snplist$CHROM, snplist$POS)

#subset snplist to only include snps in data
library(dplyr)
snplist_subset <- snplist %>%
  filter(snp %in% data$snp)

#merge this subsetted list to data, so now data has the chromosome info for each scaffold
data_updated <- cbind(data, snplist_subset)

#clean up a little, remove duplicate columns
data_updated <- select(data_updated, -10, -11, -12, -13) #para
#data_updated <- select(data_updated, -17, -18, -19, -20) #lmm

#save file as backup
write.csv(data_updated, paste(id,".orderchr.csv",sep=""))


#### rearranging by chromosome order

summary(data_updated$orientation)

#if data is autosomes only
data_updated$new_chr <- as.numeric(as.character(factor(data_updated$orientation, levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '33'), labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 33))))
data_updated$new_chr[is.na(data_updated$new_chr)] <- 34 #temporary label for unlocalized snps

#if data includes sex chromosomes W and Z
#data_updated$new_chr <- as.numeric(as.character(factor(data_updated$orientation, levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '33', 'W', 'Z'), labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 33, 34, 35))))
#data$new_chr[is.na(data$new_chr)] <- 36 #temporary label for unlocalized snps

#order data by chromosome, then by id
data_ordered <- data_updated[order(data_updated[,13], data_updated[,11]),] #para
#data_ordered <- data_updated[order(data_updated[,20], data_updated[,18]),] #lmm

#separate each chr, make chromosomal position, then put back together
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
chr33 <- subset(data_ordered, new_chr==33)
chr33$order_number <- 1:nrow(chr33)
chr34 <- subset(data_ordered, new_chr==34)
chr34$order_number <- 1:nrow(chr34)
#chr35 <- subset(data_ordered, new_chr==35)
#chr35$order_number <- 1:nrow(chr35)
#chr36 <- subset(data_ordered, new_chr==36)
#chr36$order_number <- 1:nrow(chr36)

chr_combined <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chr26, chr27, chr28, chr33, chr34) #autosomes
#chr_combined <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chr26, chr27, chr28, chr33, chr34, chr35, chr36) #all chromosomes

write.csv(chr_combined, paste(id,".orderchr.bpordered.csv",sep=""))


#### making manhattan plots

library(qqman)

#bslmm manhattan plot
upperPIP=0.03
bslmmtitle="spring rank (residuals), autosomes, 87 individuals"

png(paste(id,"_PIP_manhattan_4000x1200p40_ORDERED.png",sep=""), w=4000, h=1200, pointsize=40)
manhattan(chr_combined,main=bslmmtitle, chr="new_chr", bp="order_number", p="gamma", snp="CHR", logp=FALSE, cex.axis=0.9, col = c("darkorchid4", "gray65"), chrlabs=c(1:28, 33, "Un"), suggestiveline = FALSE, genomewideline = FALSE, ylim=c(0,upperPIP), ylab="PIP")
dev.off()


#lmmm manhattan plot
#upperlogp=9
#lmmtitle="spring rank, 5 covariates, autosomes, 87 individuals"

#png(paste(id,"_lmm_manhattan_4000x1200p40_ORDERED_pwald.png",sep=""), w=4000, h=1200, pointsize=40)
#manhattan(chr_combined,main=lmmtitle, chr="new_chr", bp="order_number", p="p_wald", snp="chr", cex.axis=0.9, col = c("darkorchid4", "gray65"), chrlabs=c(1:28, 33,"Un"), suggestiveline = FALSE, genomewideline = FALSE, ylim=c(0,upperlogp))
#dev.off()
