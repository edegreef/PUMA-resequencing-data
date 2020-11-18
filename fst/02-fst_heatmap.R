#!/usr/bin/Rscript
setwd("E:/Evelien's Dropbox/Dropbox/PUMA/Bioinformatics/reseq/FST")

fst <- read.csv("fst_67_quickview.csv")
fst_table <- fst[,c(2,3,4,5)]

library(pheatmap)
fst_matrix <- data.matrix(fst_table)
heatmap_test <- heatmap(fst_matrix)

colorset <- c("beige", "purple")
heatmap <- pheatmap(fst_matrix, 
                    color=colorRampPalette(colorset)(50),
                    cluster_rows=F,
                    cluster_cols=F,
                    border_color="gray", 
                    cellwidth=80, 
                    cellheight=80,
                    legend=TRUE,
                    main="FSTs with 4 pops",
                    na_col="white",
                    filename="fst_heatmap.jpeg")
