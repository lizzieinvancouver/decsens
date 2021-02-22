## Started 10 February 2021 ##
## By Lizzie, with help from Cat ##

## See also projects/misc/pep725/pep725check_decsens.R ##

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/projects/treegarden/decsens/analyses") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("hereliesboomboom")

d <- read.csv("pep_analyses/output/latitude_betpen_dailytemps.csv")
dm <- read.csv("pep_analyses/output/latitude_betpen_dailyandmst.csv")
pdat <- read.csv("pep_analyses/input/sitespeplatBBCH10phendat.csv")
sdat <- read.csv("pep_analyses/input/sitespeplatBBCH10.csv")

# Let's do some reformmating! Come on everybody!
# First, get a lat/long to pep_id lookup table
sdat$lat.lon <- paste(sdat$latnum, sdat$lonnum)
lookuppepid <- subset(sdat, select=c("PEP_ID", "lat.lon"))                   
pdat <- merge(pdat, lookuppepid, by="PEP_ID")

# Next, clean up dates ....
d$date <- as.Date(d$Date, format=("%Y-%m-%d"))
d$doy <- as.numeric(format(d$date, "%j"))
pdat$date <- as.Date(paste(pdat$YEAR, pdat$DAY, sep="-"), format=("%Y-%j"))
pdat$doy <- as.numeric(format(pdat$date, "%j"))

# Now do some climate summaries ...
hist(as.numeric(pdat$doy)) # hmm, Apr-May is probably best for set window, but still misses some so do Apr-June

dsummdoy <- aggregate(d["tmean"], d[c("lat", "lon", "lat.lon", "doy")], FUN=mean)

format(as.Date("1990-04-01", format=("%Y-%m-%d")), "%j")
dspcwindowall <- subset(d, doy>90 & doy<181)
dspcwindowmean <- aggregate(d["tmean"], d[c("lat", "lon", "lat.lon", "year")], FUN=mean)

windowlength <- 61
pdathere <- subset(pdat, doy>(windowlength+1))
nrow(pdathere)/nrow(pdat) # 1 row

d60d <- data.frame(lat=numeric(),
                   lon=numeric(),
                   lat.lon=character(),
                   year=numeric(),
                   doy=numeric(),
                   tmean60=numeric(),
                   n=numeric())

for (site in unique(pdat$lat.lon)){
    psub <- pdathere[which(pdathere$lat.lon==site),]
    dsub <- d[which(d$lat.lon==site),]
    for (year in unique(psub$YEAR)){
        psubyr <- psub[which(psub$YEAR==year),]
        dsubyr <- dsub[which(dsub$year==year),]
        if (nrow(psubyr)>1) {
        lastday <- mean(as.numeric(psubyr$DAY))
        } else
        lastday <- as.numeric(psubyr$DAY)
        clim60 <- dsubyr[which(dsubyr$doy<lastday & dsubyr$doy>(lastday-windowlength)),]
        d60dadd <- data.frame(
                   lat=dsubyr$lat[1],
                   lon=dsubyr$lon[1],
                   lat.lon=site,
                   year=year,
                   doy=lastday,
                   tmean60=mean(clim60$tmean, na.rm=TRUE),
                   n=nrow(psubyr))
        d60d <- rbind(d60d, d60dadd)
    }
}

dprepforstat <- merge(d60d, dspcwindowmean, by=c("lat", "lon", "lat.lon", "year"))
dim(dprepforstat)
dim(d60d) # good
dim(dspcwindowmean) # fine, as we more years of climate than phenology data

plot(dprepforstat$tmean60, dprepforstat$tmean, xlab="T mean 60 d before event", ylab="T mean static window")

# Need to add 10 if we want to use log ... so add 10 everywhere
dprepforstat$tmeanplus10 <- dprepforstat$tmean+10
dprepforstat$tmean60plus10 <- dprepforstat$tmean60+10

statz <- data.frame(lat=numeric(),
                   lon=numeric(),
                   lat.lon=character(),
                   lmslopeSpcWindow=numeric(),
                   logslopeSpcWindow=numeric(),
                   lmslope60d=numeric(),
                   logslope60d=numeric(),
                   n=numeric())

for (site in unique(pdat$lat.lon)){
    subby <- dprepforstat[which(dprepforstat$lat.lon==site),]
    statzadd <- data.frame(
                   lat=subby$lat[1],
                   lon=subby$lon[1],
                   lat.lon=site,
                   lmslopeSpcWindow=coef(lm(doy~tmeanplus10, data=subby))[2],
                   logslopeSpcWindow=coef(lm(log(doy)~log(tmeanplus10), data=subby))[2],
                   lmslope60d =coef(lm(doy~tmean60plus10, data=subby))[2],
                   logslope60d=coef(lm(log(doy)~log(tmean60plus10), data=subby))[2],
                   n=nrow(subby))
    statz <- rbind(statz, statzadd)
 }


pchhere <- 16
colz <- c("dodgerblue", "darkred", "darkorchid1")

par(mfrow=c(2,2))
plot(lmslope60d~lat, data=statz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean 60 d before event")
abline(lm(lmslope60d~lat, data=statz), col=colz[1])
summary(lm(lmslope60d~lat, data=statz))

plot(logslope60d~lat, data=statz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope w log: T mean 60 d before event")
abline(lm(logslope60d~lat, data=statz), col=colz[2])
summary(lm(logslope60d~lat, data=statz))

plot(lmslopeSpcWindow~lat, data=statz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean of Apr-June")
abline(lm(lmslopeSpcWindow~lat, data=statz), col=colz[1])
summary(lm(lmslopeSpcWindow~lat, data=statz))

plot(logslopeSpcWindow~lat, data=statz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean of Apr-June")
abline(lm(logslopeSpcWindow~lat, data=statz), col=colz[2])
summary(lm(logslopeSpcWindow~lat, data=statz))

## try to look at on similar scale ...
par(mfrow=c(1,2))
plot(lmslope60d~lat, data=statz, col=colz[1], pch=pchhere)
abline(lm(lmslope60d~lat, data=statz), col=colz[1])
summary(lm(lmslope60d~lat, data=statz))

points((logslope60d*10)~lat, data=statz, col=colz[2], pch=pchhere)
abline(lm((logslope60d*10)~lat, data=statz), col=colz[2])

# And ... are the two metrics more similar in log scale?
# Seems it ... but scaling is an issue again
plot(logslope60d~lat, data=statz, ylim=c(-0.8, 0.8))
abline(lm(logslope60d~lat, data=statz), col=colz[2])
abline(lm(logslopeSpcWindow~lat, data=statz), col=colz[3])

## Some quick plotting ...
library(ggplot2)

ggplot(dsummdoy, aes(x=as.numeric(doy), y=tmean)) +
    geom_line(color="dodgerblue") +
    facet_wrap(.~as.factor(lat.lon)) +
    xlab("Day of year") +
    ylab(expression(paste("Daily mean temperature (", degree, "C)"), sep="")) +
    theme_minimal()

ggplot(dsummdoy, aes(x=as.numeric(doy), y=tmean, group=lat, color=lat)) +
    geom_line() +
    xlab("Day of year") +
    ylab(expression(paste("Daily mean temperature (", degree, "C)"), sep="")) +
    theme_minimal()


ggplot(d, aes(x=as.numeric(doy), y=tmean, color=lat)) +
    geom_line() +
    facet_wrap(.~year) +
    xlab("Day of year") +
    ylab(expression(paste("Daily mean temperature (", degree, "C)"), sep="")) +
    theme_minimal()


## Now try some simulations using the real climate data ...
fstar <- 200

simleafout <- data.frame(lat=numeric(),
                   lon=numeric(),
                   lat.lon=character(),
                   year=numeric(),
                   leafoutdoy=numeric(),
                   tmean60=numeric())

for(asite in unique(d$lat.lon)){
    onesitedf <- d[which(d$lat.lon==asite),]
    for(oneyear in unique(onesitedf$year)){
        oneyeardf <-  onesitedf[which(onesitedf$year==oneyear),]
        oneyeardf$tmeanforgdd <- ifelse(oneyeardf$tmean>0, oneyeardf$tmean, 0)
        leafoutdoy <-  min(which(cumsum(oneyeardf$tmeanforgdd) > fstar))
        clim60 <- oneyeardf[which(oneyeardf$doy<leafoutdoy & oneyeardf$doy>(leafoutdoy-windowlength)),]
        simloadd <- data.frame(
                   lat=oneyeardf$lat[1],
                   lon=oneyeardf$lon[1],
                   lat.lon=asite,
                   year=oneyear,
                   leafoutdoy=leafoutdoy,
                   tmean60=mean(clim60$tmean, na.rm=TRUE))
        simleafout <- rbind(simleafout, simloadd)
    }
}

simleafoutwmean <- merge(dspcwindowmean, simleafout, by=c("lat", "lon", "lat.lon", "year"))
simleafoutwmean$tmeanplus10 <- simleafoutwmean$tmean+10
simleafoutwmean$tmean60plus10 <- simleafoutwmean$tmean60+10


simstatz <- data.frame(lat=numeric(),
                   lon=numeric(),
                   lat.lon=character(),
                   lmslopeSpcWindow=numeric(),
                   logslopeSpcWindow=numeric(),
                   lmslope60d=numeric(),
                   logslope60d=numeric())

for (site in unique(simleafoutwmean$lat.lon)){
    subby <- simleafoutwmean[which(simleafoutwmean$lat.lon==site),]
    statzadd <- data.frame(
                   lat=subby$lat[1],
                   lon=subby$lon[1],
                   lat.lon=site,
                   lmslopeSpcWindow=coef(lm(leafoutdoy~tmeanplus10, data=subby))[2],
                   logslopeSpcWindow=coef(lm(log(leafoutdoy)~log(tmeanplus10), data=subby))[2],
                   lmslope60d =coef(lm(leafoutdoy~tmean60plus10, data=subby))[2],
                   logslope60d=coef(lm(log(leafoutdoy)~log(tmean60plus10), data=subby))[2],
                   n=nrow(subby))
    simstatz <- rbind(simstatz, statzadd)
 }

pchhere <- 16
colz <- c("dodgerblue", "darkred", "darkorchid1")

par(mfrow=c(2,2))
plot(lmslope60d~lat, data=simstatz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean 60 d before event")
abline(lm(lmslope60d~lat, data=simstatz), col=colz[1])
summary(lm(lmslope60d~lat, data=simstatz))

plot(logslope60d~lat, data=simstatz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope w log: T mean 60 d before event")
abline(lm(logslope60d~lat, data=simstatz), col=colz[2])
summary(lm(logslope60d~lat, data=simstatz))

plot(lmslopeSpcWindow~lat, data=simstatz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean of Apr-June")
abline(lm(lmslopeSpcWindow~lat, data=simstatz), col=colz[1])
summary(lm(lmslopeSpcWindow~lat, data=simstatz))

plot(logslopeSpcWindow~lat, data=simstatz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean of Apr-June")
abline(lm(logslopeSpcWindow~lat, data=simstatz), col=colz[2])
summary(lm(logslopeSpcWindow~lat, data=simstatz))


# Compare real and sim data ...
drealsim <- merge(dprepforstat, simleafoutwmean, by=c("lat", "lon", "lat.lon", "year"),
    all.x=TRUE, suffixes=c(".real", ".sim"))
statzrealsim <- merge(statz, simstatz, by=c("lat", "lon", "lat.lon"),
    all.x=TRUE, suffixes=c(".real", ".sim"))

plot(doy~leafoutdoy, data=drealsim) # should color code by lat
abline(0,1)
plot(tmean.real~tmean.sim, data=drealsim) # safety-check!
plot(tmean60.real~tmean60.sim, data=drealsim) # should color code by lat

plot(lmslopeSpcWindow.real~lmslopeSpcWindow.sim, data=statzrealsim)
abline(0,1)
plot(logslopeSpcWindow.real~logslopeSpcWindow.sim, data=statzrealsim) 
abline(0,1)
plot(lmslope60d.real~lmslope60d.sim, data=statzrealsim)
abline(0,1)
plot(logslope60d.real~logslope60d.sim, data=statzrealsim) 
abline(0,1)
