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

# Use simplest examples for ms?
pdf("..//..//..//..//decsens/analyses/figures/ospreeforcechill.pdf", width=10, height=6)
grid.arrange(char2raw, char2logy, char2logyx, thielraw, thiellogy, thiellogyx, nrow = 2)
dev.off()
