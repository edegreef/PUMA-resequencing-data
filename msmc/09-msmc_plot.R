# plotting msmc2 outputs

library(ggplot2)

mu <- 3e-9
gen <- 2

setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/Evelien files/genomics/resequencing data/msmc/data")
PM141 <-read.table("PM141.min1MB.msmc2.final.txt", header=TRUE)
PM150 <-read.table("PM150.min1MB.msmc2.final.txt", header=TRUE)
PM223 <-read.table("PM223.min1MB.msmc2.final.txt", header=TRUE)
PMREF <-read.table("PMREF.min1MB.msmc2.final.txt", header=TRUE)

setwd("boot")

# load boot outputs
# create list of file names to load (100 per group is too much to load individually)
PM141_filenames <- list.files(pattern="PM141*")
PM150_filenames <- list.files(pattern="PM150*")
PM223_filenames <- list.files(pattern="PM223*")
PMREF_filenames <- list.files(pattern="PMREF*")

# import them together (by group)
PM141_boots <- lapply(PM141_filenames,function(i){
  read.delim(i, header=TRUE)
})

PM150_boots <- lapply(PM150_filenames,function(i){
  read.delim(i, header=TRUE)
})

PM223_boots <- lapply(PM223_filenames,function(i){
  read.delim(i, header=TRUE)
})

PMREF_boots <- lapply(PMREF_filenames,function(i){
  read.delim(i, header=TRUE)
})


# merge files from a large list format into a dataframe
PM141_boot <- do.call(rbind, PM141_boots)
PM150_boot <- do.call(rbind, PM150_boots)
PM223_boot <- do.call(rbind, PM223_boots)
PMREF_boot <- do.call(rbind, PMREF_boots)

# every 58 rows is new boot. need column to label this as run number.
PM141_boot$run <- rep(c(1:20), each=32)
PM150_boot$run <- rep(c(1:20), each=32)
PM223_boot$run <- rep(c(1:20), each=32)
PMREF_boot$run <- rep(c(1:20), each=32)

PM141_boot$location <- "P.s. arboricola (BC)"
PM150_boot$location <- "P.s. hesperia (AZ)"
PM223_boot$location <- "P.s. arboricola (CA)"
PMREF_boot$location <- "P.s. subis (MB)"

# merge
runs <- rbind(PM141_boot, PM150_boot, PM223_boot, PMREF_boot)

# add unique id for label+run
runs$plot_id <- paste(runs$location, runs$run, sep="_")
#cols <- c("PM141"="blue", "PM150"="blue", "PM223"="green", "PMREF"="black")


# putting samples together for coloring stuff
PM141$location <- "P.s. arboricola (BC)"
PM150$location <- "P.s. hesperia (AZ)"
PM223$location <- "P.s. arboricola (CA)"
PMREF$location <- "P.s. subis (MB)"

all <- rbind(PM141, PM150, PM223, PMREF)

group_colors <- c(`P.s. arboricola (BC)`="#0072B2",`P.s. hesperia (AZ)`="#F0E442", `P.s. arboricola (CA)`="#56B4E9", `P.s. subis (MB)`="#AA4371")

#boot
boot_col <- c(PM141="#F0E442",PM151="#0072B2", PM223="#56B4E9", PMREF="#AA4371")

#the 1x samples
ggplot()+
  geom_step(data=runs, aes(x=runs$left_time_boundary/mu*gen,y=(1/runs$lambda)/mu, group=plot_id, color=location), alpha=0.3)+
  geom_step(data=all, aes(x=all$left_time_boundary/mu*gen, y=(1/all$lambda)/mu, color=location), size=1.2)+
  theme_classic()+
  scale_x_log10(limits = c(5000,10000000)) +
  ylim(0,1000000)+
  labs(x="years ago", y="effective population size", color="Legend")+
  theme(panel.border=element_blank(), 
        axis.line=element_line(),
        legend.position=c(0.85,0.85), #(0,0)=bottom left, (1,1)=top right
        legend.background=element_rect(fill="white", color="gray30"),
        legend.text=element_text(size=10),
        #legend.title=element_text(size=10),
        legend.key.size=unit(0.8,'lines'),
        legend.title=element_blank())+
  scale_color_manual(values=group_colors)+
  scale_y_continuous(breaks=c(250000,500000,750000,1000000), labels=c("250000"="250","500000"="500", "750000"="750", "1000000"="1000"), limits=c(100000,1050000))+#, expand=c(0,0))+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e4,1e5,1e6, 1e7),trans='log10',labels=c("1e+04"=expression(10^4), "1e+05"=expression(10^5), "1e+06"=expression(10^6), "1e+07"=expression(10^7)), limits=c(10000,3000000), expand=c(0,0))+
  annotation_logticks(sides = "b")+
  xlab(expression(paste("Years ago")))

ggsave("MSMC_plot.png", width=6.5,height=3.5,dpi=1000)

