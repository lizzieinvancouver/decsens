## Started 13 Jan 2020
## By Cat - based off code by Lizzie in ospree repo ospree/analyses/bb_analysis/pep_sims/comparetopepsims.R

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

#Libraries
library(dplyr)
library(tidyr)

# Setting working directory. Add in your own path in an if statement for your file structure
setwd("~/Documents/git/decsens/analyses/pep_analyses") 

# get some data
# Betula puendula data from PEP (both have has GDD from 1 Jan to leafout)
# bp has mat from March 1st to June 1st and mat.lo is 30 days before leafout (uses tg -- aka mean -- data from E-OBS)
# bpalt is similar, but calculated uses txtm -- aka min and max (and we caculate the mean ourselves from those values) -- data from E-OBS) ... we don't use this currently 
bp <- read.csv("output/betpen_allchillsandgdds_45sites_mat_tntx_forsims.csv", header=TRUE)
bppre <- read.csv("output/bp_climatedatapre.csv")
bppost <- read.csv("output/bp_climatedatapost.csv")

bpclimall <- full_join(bppre, bppost)
bpclimall <- subset(bpclimall, select=c(yday, year, temp, spatial))
bpclimall$cc <- ifelse(bpclimall$year>=1950 & bpclimall$year<=1960, "1950-1960", "2000-2010")
names(bpclimall) <- c("doy", "year", "temp", "siteslist", "cc")
bpclimall <- bpclimall[!duplicated(bpclimall),]
nonyrs <- c(1950, 1961, 2000, 2011)
bpclimall <- bpclimall[!(bpclimall$year%in%nonyrs),]
#bpclimall <- bpclimall[!(bpclimall$doy<=60),]

lodate <- subset(bp, select=c("year", "siteslist", "cc", "lo"))
lodate <- lodate[!duplicated(lodate),]

lopersite <- data.frame()
dailytemps <- data.frame()
for(i in length(sitez)){ #i=1
  for(j in 1951:2010){ #j=1951
    lopersite <- lodate[(lodate$siteslist==i & lodate$year==j),]
    lo <- as.numeric(lopersite$lo)
    dailytemps <- bpclimall[(bpclimall$doy<=lo),]
    
  }
}

bptempandbb <- full_join(dailytemps, lodate)
bptempandbb <- bptempandbb[!duplicated(bptempandbb),]
bptempandbb$tempk <- bptempandbb$temp + 273.15

#### Now calculate temperature sensitivites but with log
estprecc.log <- lm(log(lo)~log(mat), data=subset(bp, cc=="1950-1960")) 
estpostcc.log <- lm(log(lo)~log(mat), data=subset(bp, cc=="2000-2010"))

logdiffbefore.after <- coef(estprecc.log)[2]-coef(estpostcc.log)[2]
# negative means a decline in sensitivity AFTER climate change

#### Now calculate temperature sensitivites without log to compare
estprecc <- lm(lo~mat, data=subset(bp, cc=="1950-1960")) 
estpostcc <- lm(lo~mat, data=subset(bp, cc=="2000-2010"))

diffbefore.after <- coef(estprecc)[2]-coef(estpostcc)[2]
# negative means a decline in sensitivity AFTER climate change