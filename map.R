
##########################################
#### 1) SIMPLE MAP
##########################################
# Loadpackages
library(raster)
library(rworldmap)
library(ggplot2)
library(tidyverse)
library(ggspatial)

# Set map boundary (xmin, xmax, ymin, ymax)
# Let's start with the whole world
boundary <- extent(-140, -50, 5, 60)
boundary

# Get map outlines from rworldmap. For "high" resolution, need to also install package rworldxtra
map_outline <- getMap(resolution="high")

# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=boundary) %>% fortify()

# Create a simple basemap
ggplot()+
  # Plot map land (fill), and outline (colour & size)
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="grey80", colour="black", size=0.5)

site_info <-  read.csv("C:/Users/eveli/Dropbox/PUMA/Evelien files/sample infos/PUMA_sample_info_87_full.csv", header=T)


# We can make this look nicer. For a simple map, I like the gray land fill and lighter gray outline, with the the water/sea being white. Can adjust however you like though!
map <- ggplot()+  
  # Plot map land (fill), and outline (colour & size)
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="grey95", colour="black", size=0.5)+
  # Remove extra space between map & axes
  coord_quickmap(expand=F)+
  # Axes labels
  xlab("Longitude")+
  ylab("Latitude")+
  # Adjust background theme/colour
  theme(
    axis.text = element_text(colour="black", size=12),
    axis.title = element_text(colour="black", size=14),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour="black", size=0.5),
    legend.text=element_text(size=12),
    legend.title=element_blank(),
    legend.key.size=unit(0.7, "cm"),
    legend.position="none",
    # could keep the panel grid lines if you want, I am excluding them but leaving the code here for manual adjustment
    #panel.grid.minor=element_line(colour="grey90", size=0.5),
    #panel.grid.major=element_line(colour="grey90", size=0.5)
  )+
  ggsn::scalebar(data=map_outline, dist=500, dist_unit="km", height=0.01,
                 transform=TRUE, model="WGS84",
                 location="bottomleft", anchor=c(x=-138, y=7),
                 st.bottom=FALSE, st.size=2.5, st.dist=0.015)

map



# Can also add species ranges
# Info on obtaining species range map shape files from IUCN: https://damariszurell.github.io/EEC-MGC/a4_RangeMaps.html

setwd("C:/Users/eveli/Dropbox/PUMA/post-grad/manuscripts/Sci Reports - mig/revision1")
# load shape file
species <- shapefile("redlist_species_data_1b865022-4a2b-438f-93b0-9f0447955155/data_0.shp")

breed <- species[species$SEASONAL==2,]

puma_map <- map +
  geom_polygon(data=breed, aes(x=long, y=lat, group=group), fill="#b297ceff", colour=NA, size=0.5, alpha=0.8)+
  geom_point(data=site_info, aes(x=Longitude, y=Latitude), pch=21, size=3, stroke=1, colour="black", fill="yellow")

#ggsave("PUMA_map_updated.png", width=9, height=8, dpi=1000)


erie <- shapefile("hydro_p_LakeErie/hydro_p_LakeErie.shp")
huron <- shapefile("hydro_p_LakeHuron/hydro_p_LakeHuron.shp")
mich <- shapefile("hydro_p_LakeMichigan/hydro_p_LakeMichigan.shp")
ont <- shapefile("hydro_p_LakeOntario/hydro_p_LakeOntario.shp")
stclair <- shapefile("hydro_p_LakeStClair/hydro_p_LakeStClair.shp")
sup <- shapefile("hydro_p_LakeSuperior/hydro_p_LakeSuperior.shp")

#http://mchp-appserv.cpe.umanitoba.ca/teaching/GISmanual/data.shtml
mb <- shapefile("simple_lakes_and_water/simple_lakes_and_water.shp")
map +
  geom_polygon(data=breed, aes(x=long, y=lat, group=group), fill="#b297ceff", colour=NA, size=0.5, alpha=0.9)+
  geom_polygon(data=erie, aes(x=long, y=lat, group=group), fill="white", colour="black", size=0.5)+
  geom_polygon(data=huron, aes(x=long, y=lat, group=group), fill="white", colour="black", size=0.5)+
  geom_polygon(data=mich, aes(x=long, y=lat, group=group), fill="white", colour="black", size=0.5)+
  geom_polygon(data=ont, aes(x=long, y=lat, group=group), fill="white", colour="black", size=0.5)+
  geom_polygon(data=stclair, aes(x=long, y=lat, group=group), fill="white", colour="black", size=0.5)+
  geom_polygon(data=sup, aes(x=long, y=lat, group=group), fill="white", colour="black", size=0.5)+
  geom_polygon(data=mb, aes(x=long, y=lat, group=group), fill="white", colour="black", size=0.5)+
  
  geom_point(data=site_info, aes(x=Longitude, y=Latitude), pch=21, size=3, stroke=1, colour="black", fill="yellow")+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill=NA, colour="black", size=0.5)

ggsave("PUMA_map_updated.png", width=9, height=8, dpi=1000)
