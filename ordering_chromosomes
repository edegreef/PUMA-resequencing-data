#!/usr/bin/Rscript

#First merging info about chromosome number and id for each scaffold on dataset

data <- read.table("SD.allchrom.LD.raw.89_para.csv", header=T) #load data want to add chr info to, this example using spring departure BSLMM data
scafs <- read.csv("purple_scafs_ordered_by_chicken.csv") #scaffold/chromosome info
scafs$id <- 1:nrow(scafs)
scafs$purple_scaf <- gsub("Q","",scafs$purple_scaf)

#Editing data file
data$best_chick_chr <- as.character(NA)
data$best_chick_id <- as.numeric(0)
data$orientation <- as.character(NA)

#Matching the chrom and scaff names to give the data file the chromosome info. This parts takes some time to run.
for (i in 1:length(data$CHR)){
  for (j in 1:length(scafs$chr_num)){
    
    if (data$CHROM[i] == scafs$purple_scaf[j]){
      data$best_chick_chr[i] <- as.character(scafs$chr_num[j])
      data$best_chick_id[i] <- as.numeric(scafs$id[j])
      data$orientation[i] <- as.character(scafs$ori[j])
    }
  }
}

#saving file
write.csv(data,"SD.allchrom.LD.raw.89_para_orderchr.csv") 


#Next, ordering data by chromosome, and by position within chromosomes

#need to make W, Z, U temporarily numeric, can change later.
#assigning chr33=29, chrW=30, chrZ=31, Unlocalized=32

summary(data$best_chick_chr)
data$new_chr <- as.numeric(as.character(factor(data$best_chick_chr,
                                               levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '33', 'W', 'Z'), labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31))))
data$new_chr[is.na(data$new_chr)] <- 32

data_ordered <- data[order(data[,20], data[,18]),] #order by chromosome, then by id

#separating each chr, making new chromosomal position id, then putting back together
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

#putting them back together
chr_combined <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chr23, chr24, chr25, chr26, chr27, chr28, chr29, chr30, chr31, chr32)

#save datafile that now has all the information to plot snps in order by chromosome and position in chromosome
write.csv(chr_combined, "SD.allchrom.LD.raw.89_para_orderchr_chr_and_bp_ordered.csv")
