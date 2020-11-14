#!/bin/bash

## quick check of lmm 

library(ggplot2)

data<-read.table("spring.rank.resid.autosomes.87.cov5.lmm.assoc.txt",h=T)
data$id <- 1:nrow(data)

# quick plot, will make a nicer plot later
png(filename="spring.rank.resid.autosomes.87.cov5.lmm.png")
plot(data$id,-log10(data$p_wald),col=as.factor(data$chr), main="spring ranked 87 individuals")
dev.off()

# checking to see if any SNPs are a sign
data$padj <- p.adjust(data$p_wald)
plot(data$id, data$padj)
data_adjusted <- data[data$padj != 1, ] 
write.csv(data_adjusted,"specialSNP.spring.rank.resid.autosomes.87.cov5.csv")
