library(car)
library(stringr)
#library(BiocManager)
#BiocManager::install("goseq")
library(goseq)

setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/Bioinformatics/reseq/GO")

#load all PUMA genes
genes <- read.csv("PUMA.genes.cc.csv", header=F)
names(genes)=c("scaffold","prog","type","bp_start","bp_end","X", "X1", "X2", "attributes")

#ordering chr
#scafs <- read.csv("purple_scafs_ordered_by_chicken.csv") #scaffold/chromosome info
#scafs$id <- 1:nrow(scafs)
#scafs$purple_scaf <- gsub("Q","",scafs$purple_scaf)

#genes$best_chick_chr <- as.character(NA)

#Matching the chrom and scaff names to give the data file the chromosome info. This parts takes some time to run
#for (i in 1:length(genes$scaffold)){
#  for (j in 1:length(scafs$chr_num)){
#    if (genes$scaffold[i] == scafs$purple_scaf[j]){
#      genes$best_chick_chr[i] <- as.character(scafs$chr_num[j])
#    }
#  }
#}


#loading snp data
snps <- read.csv("BSLMM_snps_pip0.1.csv", header=T)
snps_unique <- subset(snps, unique=="Y")
snps_unique$id <- seq(1,length(snps_unique[,1]))

snps_SA <- subset(snps, phen=="SA")
snps_SA$id <- seq(1,length(snps_SA[,1]))


#editing gene file
genes$match <- as.numeric(NA)
genes$DE <- as.numeric(NA)

for (ii in 1:length(genes[,1])){
  for (jj in 1:length(snps_SA[,1])){
    if ((genes$scaffold[ii] == snps_SA$scaffold[jj]) &&
        ((genes$bp_start[ii] - 100000) <= snps_SA$position[jj]) && 
        ((genes$bp_end[ii] + 100000) >= snps_SA$position[jj])) {
      genes$match[ii] <- snps_SA$id[jj]
      genes$DE[ii] <- 1
      break
    }
  }
}

genes$DE=recode(genes$DE,"NA='0'; else='1'")

genes_check <- subset(genes, DE==1)
genes <- genes_check
## run GO analysis #trying with ALL genes and snps first

#trying out getgo
#library("BiocManager")
#BiocManager::install("org.Gg.eg.db")
#library("org.Gg.eg.db")
#sgenes_for_go <- c("MEF2D", "RHBG", "HACE1", "DCN1", "RPL15", "DCBLD2", "TSC1", "RRP15", "RFX8", "MAP4K4","GDPD2", "PlsC", "EPB41L3", "SLC22A4")
#sup_genomes <- supportedGenomes() #use galGal3
#sup_ids <- supportedGeneIDs()
#go=getgo(genes_for_go, 'galGal3', 'ensGene',fetch.cats=c("GO:CC","GO:BP","GO:MF"))

## isolate gene names
split=str_split_fixed(genes$attributes, ";", 2)
split[,1]=gsub("ID=","",split[,1])

## isolate go terms
split2=str_split_fixed(genes$attributes,"Ontology_term=",2)
split2[,2]=sub("\\s+$", ",", gsub('(.{10})', '\\1 ', split2[,2])) #adding a space between each GO term ID in same row
split3=str_split_fixed(split2[,2],";",2) #removing the semi colon at end

## combine and estimate gene length
genes4=cbind(split[,1],genes[,1],genes[,4:5],split3[,1],genes[,11])
names(genes4)=c("gene","scaffold","bp_start","bp_end","GO","DE")
genes4$length=genes4$bp_end-genes4$bp_start
genes4_noblanks=genes4[!(is.na(genes4$GO) | genes4$GO==""), ]
rm(list=ls(pattern="split"))

## generate gene2cat input
gene2cat=str_split(genes4_noblanks[,5]," ")
names(gene2cat)=genes4_noblanks[,1]

## run analysis
pwf=nullp(genes4_noblanks$DE,bias.data=genes4_noblanks$length)
rownames(pwf)=genes4_noblanks$gene
go_output=goseq(pwf,gene2cat=gene2cat,test.cats=c("GO:BP"),method = "Wallenius")
go_output=goseq(pwf,gene2cat=gene2cat,method = "Sampling")

go_output$padjust_underrep=p.adjust(go_output$under_represented_pvalue,method="fdr")
go_output$padjust_over=p.adjust(go_output$over_represented_pvalue,method="fdr")
write.csv(go_output,file="enriched_gos_FD.csv")
