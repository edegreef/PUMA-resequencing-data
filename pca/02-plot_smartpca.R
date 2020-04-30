library(ggplot2)

#load evec txt from smartPCA output
evec <- read.table("smartpca.run1_67_auto.evec.txt") 

#plot PCA
PUMAPCA <- ggplot(data=evec, aes(x=evec[,2],y=evec[,3], colour=factor(V12)))+
  geom_point(size=2)+
  theme_classic()+
  ggtitle("PUMA PCA")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="eigenvector1", y="eigenvector2")
print(PUMAPCA + labs(colour="Population"))


# +geom_text(aes(label=evec$V1)) #add this to ggplot for sample ID to each point
