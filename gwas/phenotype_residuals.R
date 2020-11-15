#!/bin/bash

#regressing out covariates from phenotype, and using residuals as corrected phenotype (using this for bslmms, which cannot take in a covariates file)

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


#additional notes for when regressing covariates for first bloom and green-up data

#NPN 
data <- read.csv("./firstbloom-greenness_87_chunks_NPN_only.csv", header=T)
fit=lm(NPN_uptoblood4 ~ Colony+PC1+Year_blood, data)
corrected_phenotype=resid(fit)
NPN_uptoblood4 <- as.data.frame(corrected_phenotype)
colnames(NPN_uptoblood4) <- "NPN_uptoblood4"

fit=lm(NPN_old4shifted ~ Colony+PC1+Year_blood, data)
corrected_phenotype=resid(fit)
NPN_old4shifted <- as.data.frame(corrected_phenotype)
colnames(NPN_old4shifted) <- "NPN_old4shifted"

fit=lm(NPN_uptoblood10 ~ Colony+PC1+Year_blood, data)
corrected_phenotype=resid(fit)
NPN_uptoblood10 <- as.data.frame(corrected_phenotype)
colnames(NPN_uptoblood10) <- "NPN_uptoblood10"

fit=lm(NPN_old10shifted ~ Colony+Year_blood+PC1, data)
corrected_phenotype=resid(fit)
NPN_old10shifted <- as.data.frame(corrected_phenotype)
colnames(NPN_old10shifted) <- "NPN_old10shifted"
 
fit=lm(NPN_28yr_avg ~ Colony+Year_blood+PC1, data)
corrected_phenotype=resid(fit)
NPN_28yr_avg <- as.data.frame(corrected_phenotype)
colnames(NPN_28yr_avg) <- "NPN_28yr_avg"

NPN <- cbind(NPN_uptoblood4, NPN_old4shifted, NPN_uptoblood10, NPN_old10shifted, NPN_28yr_avg)
write.csv(NPN, "NPN_bloom_residuals.csv")



#MODIS
data <- read.csv("./firstbloom-greenness_87_chunks_MODIS_only.csv", header=T)

fit=lm(MODIS_uptoblood4 ~ Colony+PC1+Year_blood, data)
corrected_phenotype=resid(fit)
MODIS_uptoblood4 <- as.data.frame(corrected_phenotype)
colnames(MODIS_uptoblood4) <- "MODIS_uptoblood4"

fit=lm(MODIS_old4shifted ~ Colony+PC1+Year_blood, data)
corrected_phenotype=resid(fit)
MODIS_old4shifted <- as.data.frame(corrected_phenotype)
colnames(MODIS_old4shifted) <- "MODIS_old4shifted"

fit=lm(MODIS_8yr_avg ~ Colony+PC1+Year_blood, data)
corrected_phenotype=resid(fit)
MODIS_8yr_avg <- as.data.frame(corrected_phenotype)
colnames(MODIS_8yr_avg) <- "MODIS_8yr_avg"

MODIS <- cbind(MODIS_uptoblood4, MODIS_old4shifted, MODIS_8yr_avg)
write.csv(MODIS, "MODIS_green_residuals.csv")