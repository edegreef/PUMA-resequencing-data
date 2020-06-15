#!/usr/bin/Rscript

#making nice manhattan plots for BSLMM and LMM

library(qqman)

#loading BSLMM dataset with chr and position ordered completed in previous R script
para <- read.csv("FA_BSLMM_PIP_QN_chr_and_bp_ordered.csv", header=T)

png("FA_QN_PIP_manhattan_4000x1200p50_ORDERED_lines.png", w=4000, h=1200, pointsize=50)
manhattan(para,main="FA PIP QN - ordered", chr="new_chr", bp="order_number", p="gamma", snp="X", logp=FALSE, ylim=c(0,1), cex.axis=0.9, col = c("darkorchid4", "gray65"), ylab="PIP", suggestiveline=0.1, genomewideline = 0.25)
dev.off()


#loading LMM dataset with chr and position ordered completed in previous R script
data <- read.csv("FA_covMF_LMM_chr_and_bp_ordered.csv", header=T)

png("FD_covMF_QN_manhattan_4000x1200p40_ORDERED_pwald_lines_ylim10.png", w=4000, h=1200, pointsize=40)
manhattan(data,main="FD allchrom covMF QN - ordered", chr="new_chr", bp="order_number", p="p_wald", snp="X", cex.axis=0.9, col = c("darkorchid4", "gray65"), chrlabs=c(1:28, 33, "w", "Z", "Un"), suggestiveline = -log10(0.1/1807829), genomewideline = -log10(0.05/1807829), ylim=c(0,10))
dev.off()
