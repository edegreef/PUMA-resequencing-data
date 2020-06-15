#!/usr/bin/Rscript

setwd("E:/Evelien's Dropbox/Dropbox/PUMA/Bioinformatics/reseq/FST")

#trying to plot distance ~ fst to show IBD
data <- read.csv("fst_distance.csv", header=TRUE) #file with colony pairs and lat long coordinates and fst values

library(tidyverse)
library(geosphere)

distance <- matrix(nrow = nrow(data),ncol = 1) #creating empty matrix with 1 col

#distm is in longitude,latitude format. output is in meters
for (row in 1:nrow(data)){
  distance[row,1] <- distm(c(data[row,3],data[row,2]), 
                           c(data[row,5],data[row,4]))[1,1]
}

#adding distances to dataframe
data$distance <- distance
data$distancekm <- as.numeric(distance / 1000) #converting to km

write.csv(data, "fst_distance_added.csv")

#quick plot
plot(data$distancekm, data$fst)

#nicer plot
library(ggplot2)
library(ggpmisc)

IBD <- ggplot(data, aes(x=distancekm,y=fst))+
  geom_point(size=5, colour="tomato2") +
  theme_classic()+
  ggtitle("IBD fst by distance")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="distance (km)", y="fst") +
  geom_text(aes(label=col),  vjust = -1, nudge_y =0) +
  geom_smooth(method = "lm", se=FALSE, colour="black", formula=y~x) +
  stat_poly_eq(formula=y~x, aes(label=paste(..eq.label.., ..rr.label.., sep="~~~")), parse=TRUE)

IBD

##### Mantel's test to compare FST with distance, latitutde, and longitude
#first inputting matrices for each variable 
dist_data <- read.csv("pop_distance_matrix.csv")
fst_data <- read.csv("pop_fst_matrix.csv")
lat_data <- read.csv("pop_latitude_matrix.csv")
long_data <- read.csv("pop_longitude_matrix.csv")

library(ade4)

dist_dist <- dist(dist_data)
fst_dist <- dist(fst_data)
lat_dist <- dist(lat_data)
long_dist <- dist(long_data)

distance_mantel <- mantel.rtest(dist_dist, fst_dist, nrepet=9999)
latitude_mantel <- mantel.rtest(lat_dist, fst_dist, nrepet=9999)
longitude_mantel <- mantel.rtest(long_dist, fst_dist, nrepet=9999)
