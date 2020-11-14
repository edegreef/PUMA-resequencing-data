#!/usr/bin/Rscript

#matching chromosome info with scaffold id

data <- read.csv("imputed.qualityfiltered.snplist.CHROM.POS.csv", header=T) #load data want to add chr info to
scafs <- read.csv("purple_scafs_ordered_by_chicken.csv") #scaffold/chromosome info
scafs$id <- 1:nrow(scafs)
scafs$purple_scaf <- gsub("Q","",scafs$purple_scaf)

#adding columns in data file to fill in 
data$best_chick_chr <- as.character(NA)
data$best_chick_id <- as.numeric(0)
data$orientation <- as.character(NA)

#looping for each scaffold. this part can take a while to run.
for (i in 1:length(data$CHROM)){
  for (j in 1:length(scafs$chr_num)){
    
    if (data$CHROM[i] == scafs$purple_scaf[j]){
      data$best_chick_chr[i] <- as.character(scafs$chr_num[j])
      data$best_chick_id[i] <- as.numeric(scafs$id[j])
      data$orientation[i] <- as.character(scafs$ori[j])
    }
  }
}

#saving file
write.csv(data,"imputed.qualityfiltered.snplist.CHROM.POS.orderchr.csv")
