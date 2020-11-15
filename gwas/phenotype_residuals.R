#!/bin/bash

#regressing out covariates from phenotype, and using residuals as corrected phenotype (for bslmms, which cannot take in a covariates file)

#example below doing this for spring and fall migratory timing phenotypes

data <- read.csv("./PUMA_sample_87_combined_phen.csv", header=T)

fit = lm(spring_rank_sum ~ Colony+Year_deployed+Sex+Age+PC1, data)
corrected_phenotype = resid(fit)
spring_resid <- as.data.frame(corrected_phenotype)
colnames(spring_resid) <- "springrank_resid"

fit = lm(fall_rank_sum ~ Colony+Year_deployed+Sex+Age+PC1, data)
corrected_phenotype = resid(fit)
fall_resid <- as.data.frame(corrected_phenotype)
colnames(fall_resid) <- "fallrank_resid"

residuals <- cbind(spring_resid, fall_resid)
write.csv(residuals, "PUMA_sample_87_residual_phenotypes.csv")