#plotting lfmm results

data <- read.csv("lfmm_77_5cov_run3.csv")
data$log <- -log10(data$x)

#matching snps with chromosome info
snp <- read.csv("FASASD_23_snporder.csv", header=T)
data_snp <- cbind(data, snp)
data_snp <- subset(data_snp, select = -X)

#just a quick check that ordering worked
#data_order <- data_snp[order(data_snp[,3]),]

#prepparing for manhattan plot
data <- data_snp
summary(data$orientation)
data$new_chr <- as.numeric(as.character(factor(data$orientation,
                                               levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '33', 'W', 'Z'), labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31))))
data$new_chr[is.na(data$new_chr)] <- 32
data_ordered <- data[order(data[,9], data[,7]),] #by chr, then id

#separating each chr, making chromosomal position, then putting back together
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

chr_combined <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chr26, chr27, chr28, chr29, chr30, chr31, chr32)

write.csv(chr_combined, "lfmm_77_cov5_run3_orderchr.bpordered.csv")

library(qqman)

#plotting
png("LFMM_77_cov5_3runs_4000x1200p40_ORDERED.png", w=4000, h=1200, pointsize=40)
manhattan(chr_combined,main="LFMM - 5 covariates - 3 repeat runs", chr="new_chr", bp="order_number", p="x", snp="chr", cex.axis=0.9, col = c("darkgreen", "gray65"), chrlabs=c(1:28, 33, "w", "Z", "Un"), suggestiveline = FALSE, genomewideline = FALSE)
dev.off()

#plotting with threshold lines
png("LFMM_77_cov5_3runs_4000x1200p40_ORDERED_lines.png", w=4000, h=1200, pointsize=40)
manhattan(chr_combined,main="LFMM - 5 covariates - 3 repeat runs", chr="new_chr", bp="order_number", p="x", snp="chr", cex.axis=0.9, col = c("darkgreen", "gray65"), chrlabs=c(1:28, 33, "w", "Z", "Un"), suggestiveline = -log10(0.1/4615223), genomewideline = -log10(0.05/4615223))
dev.off()

#zooming in on certain spot
zoom <- subset(chr_combined, new_chr >= 18)
zoom <- subset(zoom, new_chr <=22)
manhattan(zoom,main="LFMM - 5 covariates - 2 repeat runs", chr="new_chr", bp="order_number", p="x", snp="chr", cex.axis=0.9, col = c("gray65", "darkgreen"), suggestiveline = FALSE, genomewideline = FALSE)