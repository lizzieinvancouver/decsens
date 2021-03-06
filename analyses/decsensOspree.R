## Started 24 February 2021 ##
## By Lizzie ##

## Trying to find a controlled multi-force and multi-chill study to use in decsens ##
## Maybe should also consider multi-photoperiod? ##

rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
library(ggplot2)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("lizzie", getwd())>0)) { 
  setwd("~/Documents/git/treegarden/budreview/ospree/bb_analysis") 
} else if (length(grep("ailene", getwd()))>0) {setwd("/Users/aileneettinger/git/ospree/analyses/bb_analysis")
}else if(length(grep("Ignacio", getwd()))>0) { 
  setwd("~/GitHub/ospree/analyses/bb_analysis") 
} else if(length(grep("catchamberlain", getwd()))>0) { 
  setwd("~/Documents/git/ospree/analyses/bb_analysis") 
} else if(length(grep("danielbuonaiuto", getwd()))>0) { 
  setwd("~/Documents/git/ospree/analyses/bb_analysis") 
} else
 setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/bb_analysis")

# get the data from ospree repo
osp <- read.csv("..//output/ospree_clean_withchill_BB.csv", header = TRUE)
dall <- subset(osp, datasetID!="junttila12") # removing junttilla which I noted is a dormancy release study (might be others that I did not remove!)
source("..//source/commoncols.R")
d <- dall[,which(colnames(dall) %in% c(common.cols.wchill, "forcetemp_night", "photoperiod_night"))]
d$datasetstudy <- paste(d$datasetID, d$study)
d$forceday <- as.numeric(d$forcetemp)
d$forcenight <- as.numeric(d$forcetemp_night)
d$photonight <- as.numeric(d$photoperiod_night)

# Get a file to estimate cue numbers 
studfile <- read.csv("..//output/studytype_table.csv", header=TRUE)
studfile <- subset(studfile, datasetID!="junttila12")

# Which study is all experimental (no ambient) and has most forcing treatments ...
manyforce <- subset(studfile, force>4)
# caffarra11a exp 2 -- okay, 5 levels force (6-32 C) x 3 chilltemp -- BIGGEST range of  temps
# campbell75 exp 3 -- okay, 4 levels force (10-22 C) x 3 chilldays at 4 C (I think) 
# charrier11 exp 2 -- okay, 5 levels force (5-25 C) x ?? chill
# fu13 is ambient
# myking97 exp1 -- okay, 6 levels force (6-21 C) x 120 days chill
# pettersen71  exp2 -- okay, 6 levels force (15-28 C) x 4 photoperiods x ?? chill

manychill <- subset(studfile, chill>7|chilltime>7)
# cannell83 exp1-2: ambient photo and NA forcetemp 
# cook00b  exp1 and cronje03  exp1 aboth (Malus domestica)
# granhus09  exp1 -- weird interrupted chilling exp
# jones12  exp2 (Ribes nigrum)
# lamb37  exp1 (Rubus nigrum)
# sonsteby14  exp1 (Ribes nigrum)
# thielges75  exp1 -- looks good for deciduous
# worrall67 exp 3 -- 2 forcetemp/photo x 13 chill time (Picea abies)

##
## So ...
## Do caffarra11a exp 2 // charrier11 exp 2 for force and thielges75  exp1 for chill
d$chilldaysnum <- as.numeric(d$chilldays)
caff2 <- d[which((d$datasetID %in% "caffarra11a") & (d$study %in% "exp2")),]
char2 <- d[which((d$datasetID %in% "charrier11") & (d$study %in% "exp2")),]
thiel <- d[which((d$datasetID %in% "thielges75") & (d$study %in% "exp1")),]

library(ggplot2)
library(gridExtra)

## forcing ...
caff22sp <- subset(caff2, genus=="Betula" | genus=="Fagus")
caffraw <- ggplot(caff22sp, aes(x=forceday, y=response.time, color=genus)) +
    geom_point() +
    xlab("Forcing temperature") +
    ylab("Days to event") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"), legend.position = "none")

cafflogy <- ggplot(caff22sp, aes(x=forceday, y=log(response.time), color=genus)) +
    geom_point() +
    xlab("Forcing temperature") +
    ylab("log(Days to event)") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"), legend.position = "none")

cafflogyx <- ggplot(caff22sp, aes(x=log(forceday), y=log(response.time), color=genus)) +
    geom_point() +
    xlab("log(Forcing temperature)") +
    ylab("log(Days to event)") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"), legend.position = "none")


char2raw <- ggplot(char2, aes(x=forceday, y=response.time)) +
    geom_point() +
    xlab("Forcing temperature") +
    ylab("Days to event") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"))

char2logy <- ggplot(char2, aes(x=forceday, y=log(response.time))) +
    geom_point() +
    xlab("Forcing temperature") +
    ylab("log(Days to event)") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"))

char2logyx <- ggplot(char2, aes(x=log(forceday), y=log(response.time))) +
    geom_point() +
    xlab("log(Forcing temperature)") +
    ylab("log(Days to event)") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"))

pdf("..//..//..//..//decsens/analyses/figures/ospreeforcecaffchar.pdf", width=10, height=6)
grid.arrange(caffraw, cafflogy, cafflogyx, char2raw, char2logy, char2logyx, nrow = 2)
dev.off()

pdf("..//..//..//..//decsens/analyses/figures/ospreeforcecharrier11.pdf", width=10, height=3)
grid.arrange(char2raw, char2logy, char2logyx, nrow = 1)
dev.off()


## chilling ...
thielraw <- ggplot(thiel, aes(x=chilldaysnum, y=response.time)) +
    geom_point() +
    xlab("Chill days") +
    ylab("Days to event") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"))

thiellogy <- ggplot(thiel, aes(x=chilldaysnum, y=log(response.time))) +
    geom_point() +
    xlab("Chill days") +
    ylab("log(Days to event)") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"))

thiellogyx <- ggplot(thiel, aes(x=log(chilldaysnum), y=log(response.time))) +
    geom_point() +
    xlab("log(Chilldays)") +
    ylab("log(Days to event)") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"))

thielUtahraw <- ggplot(thiel, aes(x=Total_Utah_Model, y=response.time)) +
    geom_point() +
    xlab("Chilling units (Utah)") +
    ylab("Days to event") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"))

thielUtahlogy <- ggplot(thiel, aes(x=Total_Utah_Model, y=log(response.time))) +
    geom_point() +
    xlab("Chilling units (Utah)") +
    ylab("log(Days to event)") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"))

thielUtahlogyx <- ggplot(thiel, aes(x=log(Total_Utah_Model), y=log(response.time))) +
    geom_point() +
    xlab("log(Chilling units (Utah))") +
    ylab("log(Days to event)") +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5, colour = "darkgray"))

pdf("..//..//..//..//decsens/analyses/figures/ospreechillthielges75.pdf", width=10, height=6)
grid.arrange(thielraw, thiellogy, thiellogyx, thielUtahraw, thielUtahlogy, thielUtahlogyx, nrow = 2)
dev.off()


##
##

# For manuscript, add sims

expected_temp <- c(rep(5,10), rep(10,10), rep(15,10), rep(20, 10), rep(25, 10))
daily_temp <- sapply(expected_temp, function(x) rnorm(100, 0 + x, 2)) # now make some simple daily temps
leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > 200))) # set leafout date as whenever 200 GDD is reached
realized_temp <- colMeans(daily_temp) # estimate the mean temp of each simulated dataset

cex.mainhere <- 1.3

pdf("..//..//..//..//decsens/analyses/figures/ospreeforcems.pdf", width=9, height=3.25)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
par(mfrow=c(1,3))
if(FALSE){
plot(jitter(leafout_date)~expected_temp, pch=19,
     ylab="Days to event",
     xlab=expression(paste("Forcing temperature (", degree, "C)")),
     main="Sims: Linear (untransformed)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
plot(jitter(log(leafout_date))~expected_temp, pch=19,
     ylab="Days to event",
     xlab=expression(paste("Forcing temperature (", degree, "C)")),
     main="Sims: Non-linear (logged y)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
plot(jitter(log(leafout_date))~log(expected_temp), pch=19,
     ylab="Days to event",
     xlab=expression(paste("Forcing temperature (", degree, "C)")),
     main="Sims: Non-linear (logged x and y)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
}
plot(response.time~forceday, data=char2, pch=19,
     ylab="Days to budburst",
     xlab=expression(paste("Forcing temperature (", degree, "C)")),
     main="Linear (untransformed)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
plot(log(response.time)~forceday, data=char2, pch=19,
     ylab="Days to budburst",
     xlab=expression(paste("Forcing temperature (", degree, "C)")),
     main="Non-linear (logged y)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
plot(log(response.time)~log(forceday), data=char2, pch=19,
     ylab="Days to budburst",
     xlab=expression(paste("Forcing temperature (", degree, "C)")),
     main="Non-linear (logged x and y)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
dev.off()


# Add a line to plots? Ailene's idea if we want to keep the sims
char2$gdd <- char2$response.time*char2$forceday
median(char2$gdd)

expected_temp <- rep(seq(from=5, to=25, by=0.05), each=10)
daily_temp <- sapply(expected_temp, function(x) rnorm(1000, 0 + x, 0.5)) # now make some simple daily temps
leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > median(char2$gdd)))) 

par(mfrow=c(1,3))
plot(response.time~forceday, data=char2, pch=19,
     ylab="Days to budburst",
     xlab=expression(paste("Forcing temperature (", degree, "C)")),
     main="Linear (untransformed)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
points(leafout_date~expected_temp, col="salmon")
points(response.time~forceday, data=char2)

plot(log(response.time)~forceday, data=char2, pch=19,
     ylab="Days to budburst",
     xlab=expression(paste("Forcing temperature (", degree, "C)")),
     main="Linear (untransformed)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
points(log(leafout_date)~expected_temp, col="salmon")
points(log(response.time)~forceday, data=char2)

plot(log(response.time)~log(forceday), data=char2, pch=19,
     ylab="Days to budburst",
     xlab=expression(paste("Forcing temperature (", degree, "C)")),
     main="Linear (untransformed)", 
     font.main = 1, cex.main = cex.mainhere, cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.25, 0), tck=-.01)
points(log(leafout_date)~log(expected_temp), col="salmon")
points(log(response.time)~log(forceday), data=char2)


