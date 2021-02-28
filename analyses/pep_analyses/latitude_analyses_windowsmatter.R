## Started 28 February 2021 ##
## By Lizzie ##

## Similar to latitude_analyses. R ## 
## But reviews some interesting window issues I found when working on this ##

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# librarires
library(ggplot2)


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

##
## Look at linear trends
## Using Swedish data only 
dprepforstatSweden <- subset(dprepforstat, lat>58)

# positive sens (higher temps -> later leafout)
plotleafouttmean60 <- ggplot(dprepforstatSweden, aes(y=as.numeric(doy), x=tmean60)) +
    geom_point() +
    geom_smooth(method = "lm", linetype = 2, lwd=0.5, color="darkgray", se = FALSE) +
    facet_wrap(.~as.factor(lat.lon)) +
    ylab("Day of year") +
    xlab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()

# negative sens (higher temps -> earlier leafout)
plotleafouttmeanAprJune <- ggplot(dprepforstatSweden, aes(y=as.numeric(doy), x=tmean)) +
    geom_point() +
    geom_smooth(method = "lm", linetype = 2, lwd=0.5, color="darkgray", se = FALSE) +
    facet_wrap(.~as.factor(lat.lon)) +
    ylab("Day of year") +
    xlab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()

pdf("pep_analyses/figures/windowsmatter/plotleafouttmean60.pdf", width=10, height=10)
plotleafouttmean60
dev.off()

pdf("pep_analyses/figures/windowsmatter/plotleafouttmeanAprJune.pdf", width=10, height=10)
plotleafouttmeanAprJune
dev.off()

# WTF? Look at the climate trends
# Now get the climate data over the 60 days ...
d60alldat <- data.frame(lat=numeric(),
                   lon=numeric(),
                   lat.lon=character(),
                   year=numeric(),
                   doy=numeric(),
                   leafout=numeric(),
                   tmean60=numeric(),
                   tmean=numeric())

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
        d60alldatadd <- data.frame(
                   lat=clim60$lat,
                   lon=clim60$lon,
                   lat.lon=clim60$lat.lon,
                   year=clim60$year,
                   doy=clim60$doy,
                   leafout=rep(lastday, nrow(clim60)),
                   tmean60=rep(mean(clim60$tmean, na.rm=TRUE), nrow(clim60)),
                   tmean=clim60$tmean)
        d60alldat <- rbind(d60alldat, d60alldatadd)
    }
}

# Let's look at two of the extreme sites ...
onesite <- subset(d60alldat, lat==69.05)
addtmean <- aggregate(onesite["tmean"], onesite[c("lat", "lon", "lat.lon", "leafout", "year")], FUN=mean)
addtmean$tmean <- round(addtmean$tmean, 2)
addtmean$label <- paste("meantemp:", addtmean$tmean, "leafout", addtmean$leafout, sep=" ")

plot69pt05 <- ggplot(onesite, aes(x=as.numeric(doy), y=tmean)) +
    geom_line() +
    facet_wrap(.~as.factor(year)) +
    geom_text(color="dodgerblue", size=3, data=addtmean, aes(x = 110, y = 13, label = label)) +
    xlab("Day of year") +
    ylab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()

nextsite <- subset(d60alldat, lat==68.0167)
addtmeannext <- aggregate(nextsite["tmean"], nextsite[c("lat", "lon", "lat.lon", "leafout", "year")], FUN=mean)
addtmeannext$tmean <- round(addtmeannext$tmean, 2)
addtmeannext$label <- paste("meantemp:", addtmeannext$tmean, "leafout", addtmeannext$leafout, sep=" ")

plot68pt02 <- ggplot(nextsite, aes(x=as.numeric(doy), y=tmean)) +
    geom_line() +
    facet_wrap(.~as.factor(year)) +
    geom_text(color="dodgerblue", size=3, data=addtmeannext, aes(x = 110, y = 13, label = label)) +
    xlab("Day of year") +
    ylab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()


pdf("pep_analyses/figures/windowsmatter/plot60dayTmean_69pt05.pdf", width=10, height=10)
plot69pt05
dev.off()

pdf("pep_analyses/figures/windowsmatter/plot60dayTmean_68pt02.pdf", width=10, height=10)
plot68pt02
dev.off()

# The issue seems to be that earlier leafout means you push the 60 days into SUPER cold snaps
# And you get colder temps on average....
# Suggests to me that accumulation of GDD is happening more rapidly now. 

## Now that we know we have THAT weirdness ...do the statistics
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

# get a list of sites with a small n of years
hist(statz$n)
smallnsites <- subset(statz, n<20)
statzsmn <- statz[which(statz$lat.lon %in% smallnsites$lat.lon),]

pchhere <- 16
colz <- c("dodgerblue", "darkred", "darkorchid1", "skyblue")

pdf("pep_analyses/figures/windowsmatter/loglinermodels_bylatitude.pdf", width=8, height=8)
par(mfrow=c(2,2))
plot(lmslope60d~lat, data=statz, col=colz[4], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean 60 d before event",
     main="Linear -- Temp 60 days before LO")
points(lmslope60d~lat, data=statzsmn, col=colz[1], pch=pchhere)
abline(lm(lmslope60d~lat, data=statz), col=colz[1])
summary(lm(lmslope60d~lat, data=statz))

plot(logslope60d~lat, data=statz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope w log: T mean 60 d before event",
     main="Logged -- Temp 60 days before LO")
abline(lm(logslope60d~lat, data=statz), col=colz[2])
summary(lm(logslope60d~lat, data=statz))

plot(lmslopeSpcWindow~lat, data=statz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean of Apr-June",
     main="Linear -- Temp 60 Apr-June")
abline(lm(lmslopeSpcWindow~lat, data=statz), col=colz[1])
summary(lm(lmslopeSpcWindow~lat, data=statz))

plot(logslopeSpcWindow~lat, data=statz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean of Apr-June",
     main="Logged -- Temp 60 days before LO")
abline(lm(logslopeSpcWindow~lat, data=statz), col=colz[2])
summary(lm(logslopeSpcWindow~lat, data=statz))
dev.off()


## Which sites are outliers?
outliers <- subset(statz, logslope60d>0.4 | logslope60d< -0.2)
outliers <- outliers[which(outliers$lat<56),]
dprepforstat[which(dprepforstat$lat.lon %in% outliers$lat.lon),]

## Some quick plotting of climate ...

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


ggplot(d60d, aes(x=tmean60, color=year)) +
    geom_histogram() +
    facet_wrap(.~lat, scales="free") +
    xlab(expression(paste("60 d mean temperature (", degree, "C)"), sep="")) +
    theme_minimal()


##
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

pdf("pep_analyses/figures/windowsmatter/loglinermodels_bylatitude_sims.pdf", width=8, height=8)
par(mfrow=c(2,2))
plot(lmslope60d~lat, data=simstatz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean 60 d before event",
     main="Linear -- Temp 60 days before LO")
abline(lm(lmslope60d~lat, data=simstatz), col=colz[1])
summary(lm(lmslope60d~lat, data=simstatz))

plot(logslope60d~lat, data=simstatz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope w log: T mean 60 d before event",
     main="Logged -- Temp 60 days before LO")
abline(lm(logslope60d~lat, data=simstatz), col=colz[2])
summary(lm(logslope60d~lat, data=simstatz))

plot(lmslopeSpcWindow~lat, data=simstatz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean of Apr-June",
     main="Linear -- Temp Apr to June")
abline(lm(lmslopeSpcWindow~lat, data=simstatz), col=colz[1])
summary(lm(lmslopeSpcWindow~lat, data=simstatz))

plot(logslopeSpcWindow~lat, data=simstatz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope w log: T mean of Apr-June",
     main="Logged -- Temp Apr to June")
abline(lm(logslopeSpcWindow~lat, data=simstatz), col=colz[2])
summary(lm(logslopeSpcWindow~lat, data=simstatz))
dev.off()


# Compare real and sim data ...
drealsim <- merge(dprepforstat, simleafoutwmean, by=c("lat", "lon", "lat.lon", "year"),
    all.x=TRUE, suffixes=c(".real", ".sim"))
statzrealsim <- merge(statz, simstatz, by=c("lat", "lon", "lat.lon"),
    all.x=TRUE, suffixes=c(".real", ".sim"))

par(mfrow=c(1,1))
plot(tmean.real~tmean.sim, data=drealsim) # safety-check!
abline(0,1)

par(mfrow=c(2,2))
plot(lmslopeSpcWindow.real~lmslopeSpcWindow.sim, data=statzrealsim)
abline(0,1)
plot(logslopeSpcWindow.real~logslopeSpcWindow.sim, data=statzrealsim) 
abline(0,1)
plot(lmslope60d.real~lmslope60d.sim, data=statzrealsim)
abline(0,1)
plot(logslope60d.real~logslope60d.sim, data=statzrealsim) 
abline(0,1)

plotleafoutrealsimsbylat <- ggplot(drealsim, aes(x=leafoutdoy, y=doy, color=lat)) +
    geom_point() +
    xlab("Day of year: Simulated") +
    ylab(expression(paste("Day of year: real"), sep="")) +
    theme_minimal()

plotleafoutrealsimsbyyr <- ggplot(drealsim, aes(x=leafoutdoy, y=doy, color=year)) +
    geom_point() +
    facet_wrap(.~lat) +
    xlab("Day of year: Simulated") +
    ylab(expression(paste("Day of year: real"), sep="")) +
    theme_minimal()

pdf("pep_analyses/figures/windowsmatter/plotleafoutrealsimsbylat.pdf", width=8, height=8)
plotleafoutrealsimsbylat
dev.off()


## Mapping ... 
library("maptools")
data(wrld_simpl)
EuropeList <- c('Germany', 'France', 'Norway', 'Sweden', 'Finland', 'Poland', 'Lithuania', 'Austria', 'Switzerland')
my_map <- wrld_simpl[wrld_simpl$NAME %in% EuropeList, ]
plot(my_map)
points(drealsim$lon, drealsim$lat, pch=16, col="skyblue", cex=0.5)
points(smallnsites$lon, smallnsites$lat, pch=16, col="darkred", cex=0.5)
points(outliers$lon, outliers$lat, pch=16, col="darkblue", cex=0.5)

