library(LEA)

project = NULL
project = snmf("Imputed_all_snps4615223.77.recode.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new")

plot(project, col = "blue", pch = 19, cex = 1.2)
