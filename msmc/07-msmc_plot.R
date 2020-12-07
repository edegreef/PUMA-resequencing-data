#plotting msmc outputs. code is a bit messy, not complete yet

library(ggplot2)

mu <- 3e-9
gen <- 2

setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/post-grad/MSMC")
PM141 <-read.table("PM141.min1MB.msmc2.final.txt", header=TRUE)
PM150 <-read.table("PM150.min1MB.msmc2.final.txt", header=TRUE)
PM223 <-read.table("PM223.min1MB.msmc2.final.txt", header=TRUE)
PMREF <-read.table("PMREF.min1MB.msmc2.final.txt", header=TRUE)

AB2 <-read.table("2AB.phased.min1MB.indices.run2.msmc2.final.txt", header=TRUE)
AB3 <-read.table("3AB.phased.min1MB.indices.run2.msmc2.final.txt", header=TRUE)
AB4 <-read.table("4AB.phased.min1MB.indices.run2.msmc2.final.txt", header=TRUE)
AB5 <-read.table("5AB.phased.min1MB.indices.run2.msmc2.final.txt", header=TRUE)

FL2 <-read.table("2FL.phased.min1MB.indices.run2.msmc2.final.txt", header=TRUE)
FL3 <-read.table("3FL.phased.min1MB.indices.run2.msmc2.final.txt", header=TRUE)
FL4 <-read.table("4FL.phased.min1MB.indices.run2.msmc2.final.txt", header=TRUE)
FL5 <-read.table("5FL.phased.min1MB.indices.run2.msmc2.final.txt", header=TRUE)

#cols <- c("PM141"="blue", "PM150"="blue", "PM223"="green", "PMREF"="black")

#AB and FL
ggplot()+
  geom_step(data=PMREF, aes(x=PMREF$left_time_boundary/mu*gen, y=(1/PMREF$lambda)/mu), color="#BBBBBB", size=1.5)+
  geom_step(data=AB5, aes(x=AB5$left_time_boundary/mu*gen, y=(1/AB5$lambda)/mu), color="#A50026", size=1.5)+
  geom_step(data=FL5, aes(x=FL5$left_time_boundary/mu*gen, y=(1/FL5$lambda)/mu), color="#FDB366", size=1.5)+
  theme_classic()+
  scale_x_log10(limits = c(1000,10000000)) +
  ylim(0,1000000)+
  labs(x="years ago", y="effective population size", color="Legend")+
  theme(legend.position="right")

#the 1x samples
ggplot()+
  geom_step(data=PMREF, aes(x=PMREF$left_time_boundary/mu*gen, y=(1/PMREF$lambda)/mu), color="#BBBBBB", size=1.5)+
  geom_step(data=PM141, aes(x=PM141$left_time_boundary/mu*gen, y=(1/PM141$lambda)/mu), color="#0072B2", size=1.5)+
  geom_step(data=PM150, aes(x=PM150$left_time_boundary/mu*gen, y=(1/PM150$lambda)/mu), color="#F0E442", size=1.5)+
  geom_step(data=PM223, aes(x=PM223$left_time_boundary/mu*gen, y=(1/PM223$lambda)/mu), color="#56B4E9", size=1.5)+
  theme_classic()+
  scale_x_log10(limits = c(5000,10000000)) +
  ylim(0,1000000)+
  labs(x="years ago", y="effective population size", color="Legend")+
  theme(legend.position="right")

#



plot(PM150$left_time_boundary/mu*gen, (1/PM150$lambda)/mu, log="x",ylim=c(0,2100000),xlim=c(10000,3000000),type="n", xlab="Years ago", ylab="effective population size")
lines(PM141$left_time_boundary/mu*gen, (1/PM141$lambda)/mu, type="s", col="red", lwd=3)
lines(PM150$left_time_boundary/mu*gen, (1/PM150$lambda)/mu, type="s", col="blue", lwd=3)
lines(PM223$left_time_boundary/mu*gen, (1/PM223$lambda)/mu, type="s", col="forestgreen", lwd=3)
lines(PMREF$left_time_boundary/mu*gen, (1/PMREF$lambda)/mu, type="s", col="black", lwd=3)

lines(AB$left_time_boundary/mu*gen, (1/AB$lambda)/mu, type="s", col="orange")
lines(FL$left_time_boundary/mu*gen, (1/FL$lambda)/mu, type="s", col="purple")

legend("topright",legend=c("PM141-BC", "PM150-AZ", "PM223-CA", "PMREF-MB"), col=c("red", "blue", "forestgreen", "black"), lty=c(1,1), lwd=3)

