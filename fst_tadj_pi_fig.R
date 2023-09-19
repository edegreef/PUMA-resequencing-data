# load inputs
early_pi <- read.csv("data/early10.auto.5kb.windowed.pi.orderedchr.plotorder.good.csv", header=TRUE)
late_pi <- read.csv("data/late10.auto.5kb.windowed.pi.orderedchr.plotorder.good.csv", header=TRUE)
early_taj <- read.csv("data/early10.auto.5kb.Tajima.D.orderedchr.plotorder.good.csv", header=TRUE)
late_taj <- read.csv("data/late10.auto.5kb.Tajima.D.orderedchr.plotorder.good.csv", header=TRUE)


# make unique ID for CHROM+BINSTART so we can merge the datasets later while keeping position
early_pi$id <- paste(early_pi$CHROM, early_pi$BIN_START, sep="-")
late_pi$id <- paste(late_pi$CHROM, late_pi$BIN_START, sep="-")

# extract specific cols for merging later
early_pi_set <- early_pi[,c(2,3,4,5,6,14)]
colnames(early_pi_set)[4] <- "early_nvariants"
colnames(early_pi_set)[5] <- "early_PI"

late_pi_set <- late_pi[,c(5,6,14)]
colnames(late_pi_set)[1] <- "late_nvariants"
colnames(late_pi_set)[2] <- "late_PI"

merge_pi <- merge(early_pi_set, late_pi_set, by="id", all=TRUE)

scafs <- read.csv("purple_scafs_ordered_by_chicken.csv", header=TRUE)
scafs$order_id <- 1:nrow(scafs)
good_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "33", "W", "Z")

scaf_order <- scafs %>% arrange(factor(chr_num, levels=good_order))
scaf_order$order_id <- 1:nrow(scafs)


# add chr info 
multi_pi <- merge_pi %>% left_join(scaf_order, by=c("CHROM"="purple_scaf")) # add chr info

# add "chru" for unplaced scaffolds
#multi_pi$target_chr <- replace_na(multi_pi$target_chr, "chru")

# removing unlocalized scafs
multi_pi <- multi_pi[!is.na(multi_pi$chr_num), ]

# order dataframe to prep for plotting by chromosome order
# order by order_id, then by BIN_START
multi_pi <- multi_pi[order(multi_pi[,"order_id"], multi_pi[,"BIN_START"]),] 
# add overall order id
multi_pi$overall_order <- 1:nrow(multi_pi) 

####stop
multi_pi_subset <- subset(multi_pi,overall_order > 5205 & overall_order < 6612)

# for adding a legend manually need to do some stupid things
early_part <- multi_pi_subset[, c(1,6,15)]
colnames(early_part)[2] <- "PI"
early_part$phen <- "early"

late_part <- multi_pi_subset[, c(1,8,15)]
colnames(late_part)[2] <- "PI"
late_part$phen <- "late"

binded <- rbind(early_part, late_part)
#cols <- c("early"="#882255","late"="#5d4dab")

full_pi <-ggplot()+
  annotate("rect", xmin = 5546, xmax =6265, ymin=-0.0005, ymax=0.0055, alpha = 0.1, fill = "orange")+
  geom_line(data=binded, aes(x=overall_order, y=PI, color=phen), alpha=0.1, size=2)+
  geom_smooth(data=binded, aes(x=overall_order, y=PI, color=phen),method="loess", span=0.02, se=F, lwd=1.5)+
  scale_color_manual(values=c("#882255", "#5d4dab"))+
  theme_classic()+
  ylab("PI")+
  xlab("")+
  labs(color="extreme")+
  scale_y_continuous(limits = c(-0.0005,0.0055), expand = c(0, 0))
#  theme(
#  legend.position = c(.95, .95),
#  legend.justification = c("right", "top"),
#  legend.box.just = "right",
#  legend.margin = margin(6, 6, 6, 6)
#)

###########
# load inputs

early_taj <- read.csv("data/early10.auto.5kb.Tajima.D.orderedchr.plotorder.good.csv", header=TRUE)
late_taj <- read.csv("data/late10.auto.5kb.Tajima.D.orderedchr.plotorder.good.csv", header=TRUE)

early_taj$BIN_START <- early_taj$BIN_START+1
late_taj$BIN_START <- late_taj$BIN_START+1


# make unique ID for CHROM+BINSTART so we can merge the datasets later while keeping position
early_taj$id <- paste(early_taj$CHROM, early_taj$BIN_START, sep="-")
late_taj$id <- paste(late_taj$CHROM, late_taj$BIN_START, sep="-")

# extract specific cols for merging later
early_tajd_set <- early_taj[,c(2,3,4,5,13)]
colnames(early_tajd_set)[3] <- "early_nvariants"
colnames(early_tajd_set)[4] <- "early_tajd"

late_tajd_set <- late_taj[,c(4,5,13)]
colnames(late_tajd_set)[1] <- "late_nvariants"
colnames(late_tajd_set)[2] <- "late_tajd"

merge_tajd <- merge_pi %>% merge(early_tajd_set, by="id", all=TRUE) %>%
  merge(late_tajd_set, by="id", all=TRUE)

scafs <- read.csv("purple_scafs_ordered_by_chicken.csv", header=TRUE)
scafs$order_id <- 1:nrow(scafs)
good_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "33", "W", "Z")

scaf_order <- scafs %>% arrange(factor(chr_num, levels=good_order))
scaf_order$order_id <- 1:nrow(scafs)


# add chr info 
multi_tajd <- merge_tajd %>% left_join(scaf_order, by=c("CHROM.x"="purple_scaf")) # add chr info

# add "chru" for unplaced scaffolds
#multi_pi$target_chr <- replace_na(multi_pi$target_chr, "chru")

# removing unlocalized scafs
multi_tajd <- multi_tajd[!is.na(multi_tajd$chr_num), ]

# order dataframe to prep for plotting by chromosome order
# order by order_id, then by BIN_START
multi_tajd <- multi_tajd[order(multi_tajd[,"order_id"], multi_tajd[,"BIN_START.x"]),] 
# add overall order id
multi_tajd$overall_order <- 1:nrow(multi_tajd) 


#multi_tajd_subset <- subset(multi_tajd,overall_order > 5328 & overall_order < 6875)
multi_tajd_subset <- subset(multi_tajd,overall_order > 5205 & overall_order < 6612)

# for adding a legend manually need to do some stupid things
early_partTD <- multi_tajd_subset[, c("id","early_tajd","overall_order")]
colnames(early_partTD)[2] <- "TajimaD"
early_partTD$phen <- "early"

late_partTD <- multi_tajd_subset[, c("id","late_tajd","overall_order")]
colnames(late_partTD)[2] <- "TajimaD"
late_partTD$phen <- "late"

binded_TD <- rbind(early_partTD, late_partTD)
#cols <- c("early"="#882255","late"="#5d4dab")

full_tajd <-ggplot()+
  annotate("rect", xmin = 5546, xmax =6265, ymin=-2.5, ymax=3, alpha = 0.1, fill = "orange")+
  geom_line(data=binded_TD, aes(x=overall_order, y=TajimaD, color=phen), alpha=0.1, size=2)+
  geom_smooth(data=binded_TD, aes(x=overall_order, y=TajimaD, color=phen),method="loess", span=0.02, se=F, lwd=1.5)+
  scale_color_manual(values=c("#882255", "#5d4dab"))+
  theme_classic()+
  ylab("Tajima's D")+
  xlab("")+
  labs(color="extreme")+
  scale_y_continuous(limits = c(-2.5,3), expand = c(0, 0))
 # theme(
#    legend.position = c(.95, .95),
##    legend.justification = c("right", "top"),
#    legend.box.just = "right",
#    legend.margin = margin(6, 6, 6, 6)
#  )



full_tajd

full_pi / full_tajd

full_tajd + theme(
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)



fst_main / full_pi / full_tajd


gene_pos <- read.csv("gene_combinedposition_yes.csv", header=T)
gene_pos_forplot <- multi_tajd_subset %>% 
  left_join(gene_pos, by=c("overall_order"="plot_order_combined"))

gene_pos_plot <- ggplot(gene_pos_forplot, aes(x=overall_order, y=gene_present))+
  geom_point(color="black", size=3, pch=19)+
  theme_minimal()+
  #ggtitle("Gene present in 5kb window")+
 theme(legend.position="none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank(),
        axis.title.x=element_blank(),
       axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank())+
  ylab("")+
  scale_y_continuous(expand = c(0, 0))

gene_pos_plot

fst_main / full_pi / full_tajd / gene_pos_plot + 
  plot_layout(heights = c(3,2,2,1))
ggsave("fst_pi_tajd_genepos_combined.png", width=12, height=9, dpi=2000)



###################### fst zoom
###########
# load inputs
AB <- read.table("data/within_AB_extremes8_5kb_notLDprune.windowed.weir.fst",header=TRUE)
FL <- read.table("data/within_FL_extremes6_5kb_notLDprune.windowed.weir.fst",header=TRUE)
main <- read.csv("data/extreme10_ABFL_fstdifference.forplot.csv",header=TRUE)


# make unique ID for CHROM+BINSTART so we can merge the datasets later while keeping position
AB$id <- paste(AB$CHROM, AB$BIN_START, sep="-")
FL$id <- paste(FL$CHROM, FL$BIN_START, sep="-")
main$id <- paste(main$CHROM, main$BIN_START, sep="-")

# extract specific cols for merging later
AB_set <- AB[,c(1,2,4,5,7)]
colnames(AB_set)[3] <- "AB_nvariants"
colnames(AB_set)[4] <- "AB_fst"

FL_set <- FL[,c(4,5,7)]
colnames(FL_set)[1] <- "FL_nvariants"
colnames(FL_set)[2] <- "FL_fst"

#main_set <- main[c(6,18,19)]
#colnames(main_set)[1] <- "main_nvariants"

#merge_all <- merge(main, AB_set, by="id", all=TRUE) %>%
#  merge(FL_set, by="id", all=TRUE)

merge_all <- main %>% left_join(AB_set, by="id") %>%
  left_join(FL_set, by="id")

#scafs <- read.csv("purple_scafs_ordered_by_chicken.csv", header=TRUE)
#scafs$order_id <- 1:nrow(scafs)
#good_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "33", "W", "Z")

#scaf_order <- scafs %>% arrange(factor(chr_num, levels=good_order))
#scaf_order$order_id <- 1:nrow(scafs)


# add chr info 
#ABFL_chr_add <- merge_all %>% left_join(scaf_order, by=c("CHROM"="purple_scaf")) # add chr info

# add "chru" for unplaced scaffolds
#multi_pi$target_chr <- replace_na(multi_pi$target_chr, "chru")

# removing unlocalized scafs
#ABFL_chr_add <- ABFL_chr_add[!is.na(ABFL_chr_add$chr_num), ]

# order dataframe to prep for plotting by chromosome order
# order by order_id, then by BIN_START
#merge_all <- merge_all[order(merge_all[,"best_chick_chr"], merge_all[,"BIN_START.x"]),] 
# add overall order id
#ABFL_chr_add$overall_order <- 1:nrow(ABFL_chr_add) 

#ABFL_chr_add_auto <- subset(ABFL_chr_add, chr_num != "Z") %>% subset(chr_num != "W")
  
#test main
ggplot(merge_all, aes(x=plot_order, y=fst_diff))+
  geom_point(aes(color=as.factor(new_chr)), alpha=0.7, size=2)+
 scale_color_manual(values=rep(c("#756bb1", "#bcbddc"), 29))+
  theme_bw()+
  ylim(c(0, 0.65))



#multi_tajd_subset <- subset(multi_tajd,overall_order > 5328 & overall_order < 6875)
merge_main_subset <- subset(merge_all,plot_order >= 5225 & plot_order <= 6548)

# for adding a legend manually need to do some stupid things
main_subset <- merge_main_subset[, c("id","fst_diff","plot_order")]
colnames(main_subset)[2] <- "fst"
main_subset$phen <- "main"

AB_fst_subset <- merge_main_subset[, c("id","AB_fst","plot_order")]
colnames(AB_fst_subset)[2] <- "fst"
AB_fst_subset$phen <- "AB"

FL_fst_subset <- merge_main_subset[, c("id","FL_fst","plot_order")]
colnames(FL_fst_subset)[2] <- "fst"
FL_fst_subset$phen <- "FL"

binded_fst <- rbind(main_subset, AB_fst_subset, FL_fst_subset)
#cols <- c("early"="#882255","late"="#5d4dab")

full_fst <-ggplot()+
 annotate("rect", xmin = 5550, xmax =6300, ymin=0, ymax=0.7, alpha = 0.1, fill = "orange")+
  geom_point(data=binded_fst, aes(x=plot_order, y=fst, color=phen), alpha=0.1, size=4)+
 #geom_smooth(data=binded_fst, aes(x=plot_order, y=fst, color=phen),method="loess", span=0.01, se=F, lwd=1.5)+
 # geom_point(data=binded_fst, aes(x=plot_order, y=fst, color=phen), alpha=0.3, size=2)+
  scale_color_manual(values=c("purple","blue", "red"))+
  theme_classic()+
  ylab("FST")+
  xlab("")+
  labs(color="group")+
  scale_y_continuous(limits = c(0,0.7), expand = c(0, 0))+
  geom_hline(yintercept=0.4, color="gray")
full_fst

full_fst <-ggplot(data=merge_main_subset)+
  annotate("rect", xmin = 5538, xmax =6235, ymin=0, ymax=0.7, alpha = 0.1, fill = "orange")+
# annotate("rect", xmin = 5550, xmax =6200, ymin=0, ymax=0.7, alpha = 0.1, fill = "orange")+
  geom_point(aes(x=plot_order, y=FL_fst), color="#882255", alpha=0.8, size=2)+
  geom_point(aes(x=plot_order, y=AB_fst), color="#8d82c4", alpha=0.8, size=2)+
 # geom_point(aes(x=plot_order, y=fst_diff), color="black", alpha=1, size=2)+

  #geom_smooth(data=binded_fst, aes(x=plot_order, y=fst, color=phen),method="loess", span=0.01, se=F, lwd=1.5)+
  # geom_point(data=binded_fst, aes(x=plot_order, y=fst, color=phen), alpha=0.3, size=2)+
 # scale_color_manual(values=c("purple","blue", "red"))+
  # "#882255", "#5d4dab"
  theme_classic()+
 # theme( legend.position="none",
         # panel.border = element_blank(),
  #       panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(), axis.text.x=element_blank(), axis.ticks = element_blank())+
  ylab("FST")+
  xlab("")+
  labs(color="group")+
  scale_y_continuous(limits = c(0,0.7), expand = c(0, 0))+
  scale_x_continuous(limits=c(5200,6600))
  #geom_hline(yintercept=0.4, color="gray")
full_fst

#ggsave("within_AB_FL_overlap_Fst_dark.png", width=8, height=3, dpi=500)
ggsave("within_AB_FL_overlap_Fst_updated_colour.png", width=11, height=3, dpi=1000)


# theme(
#    legend.position = c(.95, .95),
##    legend.justification = c("right", "top"),
#    legend.box.just = "right",
#    legend.margin = margin(6, 6, 6, 6)
#  )

