# Following the PCAdapt vignette: https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html

library(pcadapt)

setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/Bioinformatics/reseq/PCA")
snp_data <- read.pcadapt("PUMApops_auto.bed", type = "bed")

x <- pcadapt(input = snp_data, K = 20) 

#Plotting variation. Optimal K is point left of "straight" line
plot(x, option = "screeplot")

#Add population info 
#67 - there's probably a more elegant way to do this but works for my set
poplist.names <- c(rep("AB", 6),rep("PA", 1), rep("AB",1),rep("PA", 1),rep("FLD", 3), rep("PA", 1),rep("FLD", 1), rep("AB",4),rep("FLD", 1),rep("TX", 1),rep("FLD", 2),rep("PA", 1), rep("FLD",1),rep("TX", 1),rep("FLD", 2), rep("PA", 3),rep("AB", 1), rep("PA",1),rep("AB", 2),rep("PA", 1), rep("AB", 1),rep("TX", 1), rep("PA",2),rep("TX", 9),rep("AB", 3), rep("PA", 1),rep("FLD", 7), rep("PA",2),rep("FLD", 1),rep("AB", 1), rep("FLD", 1), rep("AB", 3))
print(poplist.names)

#plot PCA
plot(x, option = "scores", pop = poplist.names)

#plotting other eigenvectors
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)


########## looking at some additional things ############
#look at summary stuff, contains: scores, singular.values, loadings, zscores, af, maf, chi2.stat, gif, pvalues, pass, stat)
summary(x)
#can access elements using $ (i.e. x$pvalues)

#p-values
pvaluesdat <- as.data.frame(x$pvalues)
colnames(pvaluesdat) <- c("p")
pvaluesdat$logp <- -log10(pvaluesdat$p) 

#load SNP ID info to find out loadings by snp
snplist <- read.table("PUMApops_auto_snplist")
pvaluesdat$snp <- snplist$V1
sorted_test <- pvaluesdat[order(pvaluesdat$logp),]

#loadings for PC1
loadings <- as.data.frame(x$loadings)
loadings_pc1 <- as.data.frame(loadings$V1)
loadings_pc1$snp <- snplist$V1
colnames(loadings_pc1) <- c("loading", "snp")
sorted_loadings <- loadings_pc1[order(loadings_pc1$loading),]

#loadings for PC2
loadings_pc2 <- as.data.frame(loadings$V2)
loadings_pc2$snp <- snplist$V1
colnames(loadings_pc2) <- c("loading", "snp")
sorted_loadings2 <- loadings_pc2[order(loadings_pc2$loading),]


#Plots
plot(x, option="manhattan") #displays -log10 of p-values
plot(x, option="qqplot") #check expected uniform distribution of p-values using QQ plot
hist(x$pvalues, xlab="p-values", main=NULL, breaks=50, col="orange") #histogram of p-values
plot(x, option="stat.distribution") #

#installing qvalue
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("qvalue")

library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1 #expected false discovery rate < 10%
outliers <- which(qval < alpha)
length(outliers) #332891

#Benjamini-Hochberg
padj <- p.adjust(x$pvalues, method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
