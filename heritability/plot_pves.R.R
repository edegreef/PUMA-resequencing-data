# plotting pve, pge, and n_gamma from GEMMA outputs using the hyp and para files

library(ggplot2)

# MODIS - 4 year chunk
# older years
setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/post-grad/gemma_residuals/modis_old4shifted")
gemma.files <- list.files(pattern="modis_old4shifted.87.ZWoutall")
hyp.files<-grep(".hyp.txt",gemma.files,value=TRUE)
para.files<-grep(".param.txt",gemma.files,value=TRUE)

for (i in seq_along(hyp.files)) {
  if (i == 1) {
    hyp.dat.2 <- read.table(hyp.files[i], header=TRUE)
  } else {
    dat <- read.table(hyp.files[i], header=TRUE)
    hyp.dat.2 <- rbind(hyp.dat.2,dat)
  }
}

modis_old4shifted <- hyp.dat.2
modis_old4shifted$group <- "modis_old4shifted"

# newer years
setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/post-grad/gemma_residuals/modis_yearb4")
gemma.files <- list.files(pattern="modis_yearb4.87.ZWoutall")
hyp.files<-grep(".hyp.txt",gemma.files,value=TRUE)
para.files<-grep(".param.txt",gemma.files,value=TRUE)

for (i in seq_along(hyp.files)) {
  if (i == 1) {
    hyp.dat.2 <- read.table(hyp.files[i], header=TRUE)
  } else {
    dat <- read.table(hyp.files[i], header=TRUE)
    hyp.dat.2 <- rbind(hyp.dat.2,dat)
  }
}

modis_yearb4<-hyp.dat.2
modis_yearb4$group <- "modis_yearb4"

# merge the two
modis <- rbind(modis_old4shifted, modis_yearb4)

# plot
png("modis_pve.png", width=1000, height=1000, res=220)
ggplot(data=modis,aes(x=pve, fill=group, color=group))+
  geom_density(size=1, alpha=0.3)+
  theme_bw()
dev.off()

png("modis_pge.png", width=1000, height=1000, res=220)
ggplot(data=modis,aes(x=pge, fill=group, color=group))+
  geom_density(size=1, alpha=0.3)+
  theme_bw()
dev.off()

ggplot(data=modis,aes(x=n_gamma, fill=group, color=group))+
  geom_density(size=1, alpha=0.3)+
  theme_bw()



### Same thing again but using the NPN data (4 year chunks)
# older years
setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/post-grad/gemma_residuals/npn_old4shifted")
gemma.files <- list.files(pattern="NPN_old4shifted.63.ZWoutall")
hyp.files<-grep(".hyp.txt",gemma.files,value=TRUE)
para.files<-grep(".param.txt",gemma.files,value=TRUE)

for (i in seq_along(hyp.files)) {
  if (i == 1) {
    hyp.dat.2 <- read.table(hyp.files[i], header=TRUE)
  } else {
    dat <- read.table(hyp.files[i], header=TRUE)
    hyp.dat.2 <- rbind(hyp.dat.2,dat)
  }
}

npn_old4shifted <- hyp.dat.2
npn_old4shifted$group <- "npn_old4shifted"

# newer years
setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/post-grad/gemma_residuals/npn_yearb4")
gemma.files <- list.files(pattern="NPN_yearb4.63.ZWoutall")
hyp.files<-grep(".hyp.txt",gemma.files,value=TRUE)
para.files<-grep(".param.txt",gemma.files,value=TRUE)

for (i in seq_along(hyp.files)) {
  if (i == 1) {
    hyp.dat.2 <- read.table(hyp.files[i], header=TRUE)
  } else {
    dat <- read.table(hyp.files[i], header=TRUE)
    hyp.dat.2 <- rbind(hyp.dat.2,dat)
  }
}

npn_yearb4<-hyp.dat.2
npn_yearb4$group <- "npn_yearb4"

# merge the two
npn <- rbind(npn_old4shifted, npn_yearb4)

# plot
png("npn_pve.png", width=1000, height=1000, res=220)
ggplot(data=npn,aes(x=pve, fill=group, color=group))+
  geom_density(size=1, alpha=0.3)+
  theme_bw()
dev.off()

png("npn_pge.png", width=1000, height=1000, res=220)
ggplot(data=npn,aes(x=pge, fill=group, color=group))+
  geom_density(size=1, alpha=0.3, adjust=3)+
  theme_bw()
dev.off()

ggplot(data=npn,aes(x=n_gamma, fill=group, color=group))+
  geom_density(size=1, alpha=0.3, adjust=5)+
  theme_bw()

