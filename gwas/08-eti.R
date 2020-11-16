#!/bin/bash

#calculating equal-tailed interval from bslmm outputs from gemmma

library(bayestestR)

setwd("~/fall_resid")
id="geobirds87.ZWoutall.fallresid"

gemma.files <- list.files(pattern="geobirds87.ZWoutall.fallresid")
hyp.files<-grep(".hyp.txt",gemma.files,value=TRUE)
para.files<-grep(".param.txt",gemma.files,value=TRUE)

## look at run results
for (i in seq_along(hyp.files)) {
  if (i == 1) {
    hyp.dat.2 <- read.table(hyp.files[i], header=TRUE)
  } else {
    dat <- read.table(hyp.files[i], header=TRUE)
    hyp.dat.2 <- rbind(hyp.dat.2,dat)
  }
}

#summarize pve, pge, and n_gamma
eti(hyp.dat.2$pve)
eti(hyp.dat.2$pge)
eti(hyp.dat.2$n_gamma)
