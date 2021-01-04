## Started 31 December 2020 ##
## By Lizzie ##

## Trying to analyze output of sims_sw.R ##
## TO DO ...
# (1) Add sum in addition to aggregate statistic 1-4
# (2) Check what data it has determined is missing (and maybe fill it ourselves for speediness) #
# (3) Use real climate data but we generate leafout using a thermal sum model (maybe ... just do 300-500 GDD after Jan 1)?

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(climwin)
require(lubridate)
require(dplyr)
require(tidyr)

runSW <- FALSE # set to true if you want to run the sliding windows code (which is slow)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/treegarden/decsens/analyses/pep_analyses") 
} else setwd("~/Documents/git/decsens/analyses/pep_analyses")

# read in phenology data
namesphen <-  c("spatial", "year", "gddreq", "bbdoy")
dfbb0full <- read.csv("simmonds_slidingwin/output/fakebb_nochill_0deg.csv", col.names=namesphen)
dfbb2full <- read.csv("simmonds_slidingwin/output/fakebb_nochill_2deg.csv", col.names=namesphen)
dfbb0wchillfull <- read.csv("simmonds_slidingwin/output/fakebb_wchill_0deg.csv", col.names=namesphen)
dfbb2wchillfull <- read.csv("simmonds_slidingwin/output/fakebb_wchill_2deg.csv", col.names=namesphen)

# read in climate data
namesclim <- c("spatial", "year", "doy", "daily_temp")
dfclim0full <- read.csv("simmonds_slidingwin/output/fakeclim_nochill_0deg.csv", col.names=namesclim)
dfclim2full <- read.csv("simmonds_slidingwin/output/fakeclim_nochill_2deg.csv", col.names=namesclim)
dfclim0wchillfull <- read.csv("simmonds_slidingwin/output/fakeclim_wchill_0deg.csv", col.names=namesclim)
dfclim2wchillfull <- read.csv("simmonds_slidingwin/output/fakeclim_wchill_2deg.csv", col.names=namesclim)

# do some data look-sees
range(dfbb0full$bbdoy)
range(dfbb2full$bbdoy)
range(dfbb0wchillfull$bbdoy) # latest date is 171, maybe use late-June as refday?
range(dfbb2wchillfull$bbdoy)

par(mfrow=c(2,2))
hist(dfbb0full$gddreq)
hist(dfbb2full$gddreq)
hist(dfbb0wchillfull$gddreq)
hist(dfbb2wchillfull$gddreq)

##
## Data formatting 

# Need a four digit year ... so making up a 'baseyear'
baseyear <- 1900

formatphen <- function(dffull, baseyear){
    df <- dffull[which(dffull[["year"]]>1),] # delete first bb observation since I did not build climate data for a year before it
    df[["Year"]] <- df[["year"]] + baseyear
    df[["gddreq"]] <- NULL
    df[["bb_date"]] <- as.Date(paste(df[["Year"]], df[["bbdoy"]], sep="-"), format="%Y-%j")
    df[["bb_mean"]] <- mean(df[["bbdoy"]])
    return(df)
}

formatclim <- function(dffull, baseyear){
    df <- dffull
    df[["Year"]] <- df[["year"]] + baseyear
    df[["date"]] <- as.Date(paste(df[["Year"]], df[["doy"]], sep="-"), format="%Y-%j")
    return(df)
}

dfbb0 <- formatphen(dfbb0full, baseyear)
dfbb2 <- formatphen(dfbb2full, baseyear)
dfbb0wchill <- formatphen(dfbb0wchillfull, baseyear)
dfbb2wchill <- formatphen(dfbb2wchillfull, baseyear)

dfclim0 <- formatclim(dfclim0full, baseyear)
dfclim2 <- formatclim(dfclim2full, baseyear)
dfclim0wchill <- formatclim(dfclim0wchillfull, baseyear)
dfclim2wchill <- formatclim(dfclim2wchillfull, baseyear)

if(FALSE){
par(mfrow=c(1,1))
dfhere <- dfclim2wchill
plot(dfhere$daily_temp[1:1600]~c(1:1600), type="l")
}

##
## Run the sliding windows ...

runSW_nospatial <- function(datafile, climate, refday){
    results <- slidingwin(
      baseline = lm(datafile$bb_mean ~ 1), # baseline model.
      xvar=list(climate$daily_temp), 
      cdate = climate$date, 
      bdate = datafile$bb_date,
      cmissing="method1", cinterval="day", type="absolute",
      range=c(366,0), refday = refday,
      stat=c("mean", "min", "max", "slope"),
      func="lin"
    )
    return(results)
}

runSW_nospatialnoNA <- function(datafile, climate, refday){
    results <- slidingwin(
      baseline = lm(datafile$bb_mean ~ 1), # baseline model.
      xvar=list(climate$daily_temp), 
      cdate = climate$date, 
      bdate = datafile$bb_date,
      cmissing=FALSE, cinterval="day", type="absolute",
      range=c(365,0), refday = refday,
      stat=c("mean", "min", "max", "slope"),
      func="lin"
    )
    return(results)
}

runSW_nospatial_rel<- function(datafile, climate, refday){
    results <- slidingwin(
      baseline = lm(datafile$bb_mean ~ 1), # baseline model.
      xvar=list(climate$daily_temp), 
      cdate = climate$date, 
      bdate = datafile$bb_date,
      cmissing="method1", cinterval="day", type="relative",
      range=c(366,0), refday = refday,
      stat=c("mean", "min", "max", "slope"),
      func="lin"
    )
    return(results)
}

refdayhere <- c(30, 06)
# Something very funny about slidingwin f(x) is making it hard for me to functionalize
# does it actually require those named objects?!
# or maybe it's just late.
if(runSW){
datafile <- subset(dfbb0, spatial==1)
climate <- subset(dfclim0, spatial==1)
results0deg <- runSW_nospatial(datafile, climate, refdayhere)
results0degrel <- runSW_nospatial_rel(datafile, climate, refdayhere)
# results0degtest <- runSW_nospatialnoNA(datafile, climate, refdayhere) # these are missing leapday values I suspect based on how many are missing 


datafile <- subset(dfbb2, spatial==1)
climate <- subset(dfclim2, spatial==1)
results2deg <- runSW_nospatial(datafile, climate, refdayhere)
results2degrel <- runSW_nospatial_rel(datafile, climate, refdayhere)


datafile <- subset(dfbb0wchill, spatial==1)
climate <- subset(dfclim0wchill, spatial==1)
results0degwchill <- runSW_nospatial(datafile, climate, refdayhere)

datafile <- subset(dfbb2wchill, spatial==1)
climate <- subset(dfclim2wchill, spatial==1)
results2degwchill <- runSW_nospatial(datafile, climate, refdayhere)

saveRDS(results0deg, "simmonds_slidingwin/output/results0deg.RDS")
saveRDS(results2deg, "simmonds_slidingwin/output/results2deg.RDS")

saveRDS(results0degrel, "simmonds_slidingwin/output/results0degrel.RDS")
saveRDS(results2degrel, "simmonds_slidingwin/output/results2degrel.RDS")

saveRDS(results0degwchill, "simmonds_slidingwin/output/results0degwchill.RDS")
saveRDS(results2degwchill, "simmonds_slidingwin/output/results2degwchill.RDS")
}

if(!runSW){
results0deg <- readRDS("simmonds_slidingwin/output/results0deg.RDS")
results2deg <- readRDS("simmonds_slidingwin/output/results2deg.RDS")

results0degrel <- readRDS("simmonds_slidingwin/output/results0degrel.RDS")
results2degrel <- readRDS("simmonds_slidingwin/output/results2degrel.RDS")

results0degwchill <- readRDS("simmonds_slidingwin/output/results0degwchill.RDS")
results2degwchill <- readRDS("simmonds_slidingwin/output/results2degwchill.RDS")
}


if(FALSE){
## Below is not running so let's skip space for now ...
datafile <- dfbb0
climate <- dfclim0
datafile$bb_date <- strptime(datafile$bb_date, format = "%Y-%m-%d")
climate$date <- strptime(climate$date, format="%Y-%m-%d")

    results <- slidingwin(
      baseline = lm(datafile$bb_mean ~ 1), # baseline model.
      xvar=list(climate$daily_temp), 
      cdate = climate$date, 
      bdate = datafile$bb_date,
      cmissing=FALSE, cinterval="day", type="absolute",
      spatial=list(as.factor(datafile$spatial), as.factor(climate$spatial)),
      range=c(366,0), refday = refday,
      stat=c("mean", "min", "max", "slope"),
      func="lin"
    )
}

## 
## Plotting, needs clean-up!
# Clean-up on aisle 6 (ut's just copy and pasted now)
par(mfrow=c(2,2))
refdayhere
refdoy <- as.numeric(format(as.Date("30-06-1900", format="%d-%m-%Y"), format="%j"))
library("viridis")
colz <- plasma(4, alpha = 0.5)

usr <- par('usr')

## Sim data with no chill, no warming ... 
datafile <- subset(dfbb0, spatial==1)
climate <- subset(dfclim0, spatial==1)
resultshere <- results0degrel

meanclim <- aggregate(climate["daily_temp"], climate["doy"], FUN=mean)
meanclim$doyadj <- meanclim$doy-refdoy

stat1open <- refdoy-resultshere[[1]]$Dataset$WindowOpen[1]
stat1close <- refdoy-resultshere[[1]]$Dataset$WindowClose[1]
stat2open <- refdoy-resultshere[[2]]$Dataset$WindowOpen[1]
stat2close <- refdoy-resultshere[[2]]$Dataset$WindowClose[1]
stat3open <- refdoy-resultshere[[3]]$Dataset$WindowOpen[1]
stat3close <- refdoy-resultshere[[3]]$Dataset$WindowClose[1]
stat4open <- refdoy-resultshere[[4]]$Dataset$WindowOpen[1]
stat4close <- refdoy-resultshere[[4]]$Dataset$WindowClose[1]

plot(daily_temp~doyadj, data=meanclim, type="l", main="Thermal sum only, 0 C warm", xlab="Day of year adjusted to ref day")
rect(stat1open, usr[3], stat1close, usr[4], col=colz[1])
rect(stat2open, usr[3], stat2close, usr[4], col=colz[2])
rect(stat3open, usr[3], stat3close, usr[4], col=colz[3])
rect(stat4open, usr[3], stat4close, usr[4], col=colz[4])
abline(v=(datafile$bb_mean[1]-refdoy), col="seagreen4", lwd=2, lty=2)

## Sim data with no chill, 2 C warming ... 
datafile <- subset(dfbb2, spatial==1)
climate <- subset(dfclim2, spatial==1)
resultshere <- results2deg

meanclim <- aggregate(climate["daily_temp"], climate["doy"], FUN=mean)
meanclim$doyadj <- meanclim$doy-refdoy

stat1open <- refdoy-resultshere[[1]]$Dataset$WindowOpen[1]
stat1close <- refdoy-resultshere[[1]]$Dataset$WindowClose[1]
stat2open <- refdoy-resultshere[[2]]$Dataset$WindowOpen[1]
stat2close <- refdoy-resultshere[[2]]$Dataset$WindowClose[1]
stat3open <- refdoy-resultshere[[3]]$Dataset$WindowOpen[1]
stat3close <- refdoy-resultshere[[3]]$Dataset$WindowClose[1]
stat4open <- refdoy-resultshere[[4]]$Dataset$WindowOpen[1]
stat4close <- refdoy-resultshere[[4]]$Dataset$WindowClose[1]

plot(daily_temp~doyadj, data=meanclim, type="l", main="Thermal sum only, 2 C warm", xlab="Day of year adjusted to ref day")
rect(stat1open, usr[3], stat1close, usr[4], col=colz[1])
rect(stat2open, usr[3], stat2close, usr[4], col=colz[2])
rect(stat3open, usr[3], stat3close, usr[4], col=colz[3])
rect(stat4open, usr[3], stat4close, usr[4], col=colz[4])
abline(v=(datafile$bb_mean[1]-refdoy), col="seagreen4", lwd=2, lty=2)

## Sim data with chill, no warming ... 
datafile <- subset(dfbb0wchill, spatial==1)
climate <- subset(dfclim0wchill, spatial==1)
resultshere <- results0degwchill

meanclim <- aggregate(climate["daily_temp"], climate["doy"], FUN=mean)
meanclim$doyadj <- meanclim$doy-refdoy

stat1open <- refdoy-resultshere[[1]]$Dataset$WindowOpen[1]
stat1close <- refdoy-resultshere[[1]]$Dataset$WindowClose[1]
stat2open <- refdoy-resultshere[[2]]$Dataset$WindowOpen[1]
stat2close <- refdoy-resultshere[[2]]$Dataset$WindowClose[1]
stat3open <- refdoy-resultshere[[3]]$Dataset$WindowOpen[1]
stat3close <- refdoy-resultshere[[3]]$Dataset$WindowClose[1]
stat4open <- refdoy-resultshere[[4]]$Dataset$WindowOpen[1]
stat4close <- refdoy-resultshere[[4]]$Dataset$WindowClose[1]

plot(daily_temp~doyadj, data=meanclim, type="l", main="Thermal sum with chill, 0 C warm", xlab="Day of year adjusted to ref day")
rect(stat1open, usr[3], stat1close, usr[4], col=colz[1])
rect(stat2open, usr[3], stat2close, usr[4], col=colz[2])
rect(stat3open, usr[3], stat3close, usr[4], col=colz[3])
rect(stat4open, usr[3], stat4close, usr[4], col=colz[4])
abline(v=(datafile$bb_mean[1]-refdoy), col="seagreen4", lwd=2, lty=2)

## Sim data with chill, 2 C warming ... 
datafile <- subset(dfbb2wchill, spatial==1)
climate <- subset(dfclim2wchill, spatial==1)
resultshere <- results2degwchill

meanclim <- aggregate(climate["daily_temp"], climate["doy"], FUN=mean)
meanclim$doyadj <- meanclim$doy-refdoy

stat1open <- refdoy-resultshere[[1]]$Dataset$WindowOpen[1]
stat1close <- refdoy-resultshere[[1]]$Dataset$WindowClose[1]
stat2open <- refdoy-resultshere[[2]]$Dataset$WindowOpen[1]
stat2close <- refdoy-resultshere[[2]]$Dataset$WindowClose[1]
stat3open <- refdoy-resultshere[[3]]$Dataset$WindowOpen[1]
stat3close <- refdoy-resultshere[[3]]$Dataset$WindowClose[1]
stat4open <- refdoy-resultshere[[4]]$Dataset$WindowOpen[1]
stat4close <- refdoy-resultshere[[4]]$Dataset$WindowClose[1]

plot(daily_temp~doyadj, data=meanclim, type="l", main="Thermal sum with chill, 2 C warm", xlab="Day of year adjusted to ref day")
rect(stat1open, usr[3], stat1close, usr[4], col=colz[1])
rect(stat2open, usr[3], stat2close, usr[4], col=colz[2])
rect(stat3open, usr[3], stat3close, usr[4], col=colz[3])
rect(stat4open, usr[3], stat4close, usr[4], col=colz[4])
abline(v=(datafile$bb_mean[1]-refdoy), col="seagreen4", lwd=2, lty=2)


##
## Couple of plots of relative results ...

par(mfrow=c(1,2))

## Sim data with no chill, no warming ... 
datafile <- subset(dfbb0, spatial==1)
climate <- subset(dfclim0, spatial==1)
resultshere <- results0degrel

meanclim <- aggregate(climate["daily_temp"], climate["doy"], FUN=mean)
meanclim$doyadj <- meanclim$doy-refdoy

stat1open <- refdoy-resultshere[[1]]$Dataset$WindowOpen[1]
stat1close <- refdoy-resultshere[[1]]$Dataset$WindowClose[1]
stat2open <- refdoy-resultshere[[2]]$Dataset$WindowOpen[1]
stat2close <- refdoy-resultshere[[2]]$Dataset$WindowClose[1]
stat3open <- refdoy-resultshere[[3]]$Dataset$WindowOpen[1]
stat3close <- refdoy-resultshere[[3]]$Dataset$WindowClose[1]
stat4open <- refdoy-resultshere[[4]]$Dataset$WindowOpen[1]
stat4close <- refdoy-resultshere[[4]]$Dataset$WindowClose[1]

plot(daily_temp~doyadj, data=meanclim, type="l", main="Thermal sum only, 0 C warm", xlab="Day of year adjusted to ref day")
rect(stat1open, usr[3], stat1close, usr[4], col=colz[1])
rect(stat2open, usr[3], stat2close, usr[4], col=colz[2])
rect(stat3open, usr[3], stat3close, usr[4], col=colz[3])
rect(stat4open, usr[3], stat4close, usr[4], col=colz[4])
abline(v=(datafile$bb_mean[1]-refdoy), col="seagreen4", lwd=2, lty=2)

## Sim data with no chill, 2 C warming ... 
datafile <- subset(dfbb2, spatial==1)
climate <- subset(dfclim2, spatial==1)
resultshere <- results2degrel

meanclim <- aggregate(climate["daily_temp"], climate["doy"], FUN=mean)
meanclim$doyadj <- meanclim$doy-refdoy

stat1open <- refdoy-resultshere[[1]]$Dataset$WindowOpen[1]
stat1close <- refdoy-resultshere[[1]]$Dataset$WindowClose[1]
stat2open <- refdoy-resultshere[[2]]$Dataset$WindowOpen[1]
stat2close <- refdoy-resultshere[[2]]$Dataset$WindowClose[1]
stat3open <- refdoy-resultshere[[3]]$Dataset$WindowOpen[1]
stat3close <- refdoy-resultshere[[3]]$Dataset$WindowClose[1]
stat4open <- refdoy-resultshere[[4]]$Dataset$WindowOpen[1]
stat4close <- refdoy-resultshere[[4]]$Dataset$WindowClose[1]

plot(daily_temp~doyadj, data=meanclim, type="l", main="Thermal sum only, 2 C warm", xlab="Day of year adjusted to ref day")
rect(stat1open, usr[3], stat1close, usr[4], col=colz[1]) # purple
rect(stat2open, usr[3], stat2close, usr[4], col=colz[2]) # pink
rect(stat3open, usr[3], stat3close, usr[4], col=colz[3]) # pink-yellow (flesh-toned)
rect(stat4open, usr[3], stat4close, usr[4], col=colz[4]) # yellow
abline(v=(datafile$bb_mean[1]-refdoy), col="seagreen4", lwd=2, lty=2)
