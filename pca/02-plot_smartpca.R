library(ggplot2)

#load evec txt from smartPCA output
evec <- read.table("smartpca.run1_67_auto.evec.txt") 

#plot PCA
PUMAPCA <- ggplot(data=evec, aes(x=evec[,2],y=evec[,3]))+
  geom_point(size=20, alpha=0.70, aes(color=V12)) +
  theme_classic()+
  ggtitle("PUMA PCA")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="PC1", y="PC2")

#setting up grids
grids <- theme(panel.grid.major = element_line(size = (1), colour="grey"))

#printing PCA with manually added color scheme
PUMAPCA + scale_color_manual(values = c("#D55E00", "#0072B2", "#009E73", "#F0E442")) + grids

# +geom_text(aes(label=evec$V1)) #add this to ggplot for sample ID to each point
