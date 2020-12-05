#Merging .zscore outputs from lfmm models, using tutotial from: http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/note.pdf

setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/Bioinformatics/reseq/lfmm/bloom77_noprune/5cov")

z.table = NULL

#number of .zscores to merge will depend on number of covariates and run repetitions
a <- read.table("Imputed_all_snps4615223.77.recode_r1_a1.1.zscore")
b <- read.table("Imputed_all_snps4615223.77.recode_r1_a2.1.zscore")
c <- read.table("Imputed_all_snps4615223.77.recode_r1_a3.1.zscore")
d <- read.table("Imputed_all_snps4615223.77.recode_r1_a4.1.zscore")
e <- read.table("Imputed_all_snps4615223.77.recode_r1_a5.1.zscore")
f <- read.table("Imputed_all_snps4615223.77.recode_r1_a6.1.zscore")
g <- read.table("Imputed_all_snps4615223.77.recode_r7_a1.1.zscore")
h <- read.table("Imputed_all_snps4615223.77.recode_r7_a2.1.zscore")
i <- read.table("Imputed_all_snps4615223.77.recode_r7_a3.1.zscore")
j <- read.table("Imputed_all_snps4615223.77.recode_r7_a4.1.zscore")
k <- read.table("Imputed_all_snps4615223.77.recode_r7_a5.1.zscore")
l <- read.table("Imputed_all_snps4615223.77.recode_r7_a6.1.zscore")
m <- read.table("Imputed_all_snps4615223.77.recode_r13_a1.1.zscore")
n <- read.table("Imputed_all_snps4615223.77.recode_r13_a2.1.zscore")
o <- read.table("Imputed_all_snps4615223.77.recode_r13_a3.1.zscore")
p <- read.table("Imputed_all_snps4615223.77.recode_r13_a4.1.zscore")
q <- read.table("Imputed_all_snps4615223.77.recode_r13_a5.1.zscore")
r <- read.table("Imputed_all_snps4615223.77.recode_r13_a6.1.zscore")

z.table = cbind(z.table, a[,1], b[,1], c[,1], d[,1], e[,1], f[,1], g[,1], h[,1], i[,1], j[,1], k[,1], l[,1], m[,1], n[,1], o[,1], p[,1], q[,1], r[,1])

z.score = apply(z.table, MARGIN = 1, median)
lambda = median(z.score^2)/0.456
adjusted.p.values = pchisq(z.score^2/lambda, df = 1, lower = F)

hist(adjusted.p.values, col = 3)

q = 0.1
L = length(adjusted.p.values)
w = which(sort(adjusted.p.values) < q * (1:L) / L)
candidates = order(adjusted.p.values)[w]

#estimated FDR and True Positif
estimated.FDR = length(which(candidates <= 900))/length(candidates)
estimated.TP = length(which(candidates > 900))/100
print(paste("expected FDR:", q))
print(paste("FDR:", estimated.FDR, "True Positive:", estimated.TP))


write.csv(adjusted.p.values, "lfmm_77_5cov_run3.csv")
