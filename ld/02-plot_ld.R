#plotting plink LD info

library(dplyr)
library(stringr)
setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/post-grad/LD")
dfr<-read.delim("geobirds87_ld.ld",sep="",check.names=F,stringsAsFactors=F)
dfr$dist=abs(dfr$BP_B - dfr$BP_A)
dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=100))
dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(R2),median=median(R2))

plot(dfr1$mean, type="l", xlab="100 KB block", ylab="R2 mean") #line plot
plot(dfr1$mean, type="p", xlab="100 KB block", ylab="R2 mean") #dot plot
