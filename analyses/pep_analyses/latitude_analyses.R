## Started 10 February 2021 ##
## By Lizzie, with help from Cat ##

## See also projects/misc/pep725/pep725check_decsens.R ##

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
# Use only Sweden for consistent sampling ...
sdat <- subset(sdat, latnum>58)
# Get a lat/long to pep_id lookup table
sdat$lat.lon <- paste(sdat$latnum, sdat$lonnum)
lookuppepid <- subset(sdat, select=c("PEP_ID", "lat.lon"))
pdat <- merge(pdat, lookuppepid, by="PEP_ID")

# Next, clean up dates ....
d$date <- as.Date(d$Date, format=("%Y-%m-%d"))
d$doy <- as.numeric(format(d$date, "%j"))
pdat$date <- as.Date(paste(pdat$YEAR, pdat$DAY, sep="-"), format=("%Y-%j"))
pdat$doy <- as.numeric(format(pdat$date, "%j"))

# Now do some climate summaries ...
hist(as.numeric(pdat$doy)) 
range(as.numeric(pdat$doy))
format(as.Date("1990-04-15", format=("%Y-%m-%d")), "%j")
format(as.Date("1990-05-01", format=("%Y-%m-%d")), "%j")
format(as.Date("1990-05-15", format=("%Y-%m-%d")), "%j")
format(as.Date("1990-06-15", format=("%Y-%m-%d")), "%j")
# Let's do mid-April to mid-June

dsummdoy <- aggregate(d["tmean"], d[c("lat", "lon", "lat.lon", "doy")], FUN=mean)
dspcwindowall <- subset(d, doy>104 & doy<167)
dspcwindowmean <- aggregate(d["tmean"], d[c("lat", "lon", "lat.lon", "year")], FUN=mean)

windowlength <- 31
pdathere <- subset(pdat, doy>(windowlength+1))
nrow(pdathere)/nrow(pdat) 

dclimbefleaout <- data.frame(lat=numeric(),
                   lon=numeric(),
                   lat.lon=character(),
                   year=numeric(),
                   doy=numeric(),
                   tmeanLO=numeric(),
                   n=numeric())

for (site in unique(pdat$lat.lon)){
    psub <- pdathere[which(pdathere$lat.lon==site),]
    dsub <- d[which(d$lat.lon==site),]
    for (year in unique(psub$YEAR)){
        psubyr <- psub[which(psub$YEAR==year),]
        dsubyr <- dsub[which(dsub$year==year),]
        for (rowhere in 1:nrow(psubyr)){
        lastday <- as.numeric(psubyr$DAY[rowhere])
        climwindow <- dsubyr[which(dsubyr$doy<lastday & dsubyr$doy>(lastday-windowlength)),]
        dclimbefleaoutadd <- data.frame(
                   lat=dsubyr$lat[rowhere],
                   lon=dsubyr$lon[rowhere],
                   lat.lon=site,
                   year=year,
                   doy=lastday,
                   tmeanLO=mean(climwindow$tmean, na.rm=TRUE),
                   n=nrow(psubyr))
        dclimbefleaout <- rbind(dclimbefleaout, dclimbefleaoutadd)
        }
    }
}

dprepforstat <- merge(dclimbefleaout, dspcwindowmean, by=c("lat", "lon", "lat.lon", "year"))
dim(dprepforstat)
dim(dclimbefleaout) # good
dim(dspcwindowmean) # fine, as we more years of climate than phenology data

plot(dprepforstat$tmeanLO, dprepforstat$tmean, xlab="T mean 30 d before event", ylab="T mean static window")

# Need to add 10 if we want to use log ... so add 10 everywhere
dprepforstat$tmeanplus10 <- dprepforstat$tmean+10
dprepforstat$tmeanLOplus10 <- dprepforstat$tmeanLO+10

# Still weird, even with 30 day window ... 
ggplot(dprepforstat, aes(y=as.numeric(doy), x=tmeanLO)) +
    geom_point() +
    geom_smooth(method = "lm", linetype = 2, lwd=0.5, color="darkgray", se = FALSE) +
    facet_wrap(.~as.factor(lat.lon)) +
    ylab("Day of year") +
    xlab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()

ggplot(dprepforstat, aes(y=as.numeric(doy), x=tmean)) +
    geom_point() +
    geom_smooth(method = "lm", linetype = 2, lwd=0.5, color="darkgray", se = FALSE) +
    facet_wrap(.~as.factor(lat.lon)) +
    ylab("Day of year") +
    xlab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()


## Now do the statistics! 
statz <- data.frame(lat=numeric(),
                   lon=numeric(),
                   lat.lon=character(),
                   lmslopeSpcWindow=numeric(),
                   logslopeSpcWindow=numeric(),
                   lmslopeTLOd=numeric(),
                   logslopeTLOd=numeric(),
                   n=numeric())

for (site in unique(pdat$lat.lon)){
    subby <- dprepforstat[which(dprepforstat$lat.lon==site),]
    statzadd <- data.frame(
                   lat=subby$lat[1],
                   lon=subby$lon[1],
                   lat.lon=site,
                   lmslopeSpcWindow=coef(lm(doy~tmeanplus10, data=subby))[2],
                   logslopeSpcWindow=coef(lm(log(doy)~log(tmeanplus10), data=subby))[2],
                   lmslopeTLOd =coef(lm(doy~tmeanLOplus10, data=subby))[2],
                   logslopeTLOd=coef(lm(log(doy)~log(tmeanLOplus10), data=subby))[2],
                   n=nrow(subby))
    statz <- rbind(statz, statzadd)
 }

# get a list of sites with a small n of years
hist(statz$n)
smallnsites <- subset(statz, n<20)
statzsmn <- statz[which(statz$lat.lon %in% smallnsites$lat.lon),]

pchhere <- 16
colz <- c("dodgerblue", "darkred", "darkorchid1", "skyblue")

pdf("pep_analyses/figures/loglinermodels_bylatitude.pdf", width=8, height=8)
par(mfrow=c(2,2))
plot(lmslopeTLOd~lat, data=statz, col=colz[4], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean window before event",
     main="Linear -- Temp window before LO")
points(lmslopeTLOd~lat, data=statzsmn, col=colz[1], pch=pchhere)
abline(lm(lmslopeTLOd~lat, data=statz), col=colz[1])

plot(logslopeTLOd~lat, data=statz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope w log: T mean window before event",
     main="Logged -- Temp window before LO")
abline(lm(logslopeTLOd~lat, data=statz), col=colz[2])

plot(lmslopeSpcWindow~lat, data=statz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean of Apr-June",
     main="Linear -- Temp window static")
abline(lm(lmslopeSpcWindow~lat, data=statz), col=colz[1])

plot(logslopeSpcWindow~lat, data=statz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope w log: T mean of Apr-June",
     main="Logged -- Temp window static")
abline(lm(logslopeSpcWindow~lat, data=statz), col=colz[2])
dev.off()


## Some quick plotting of climate across days...

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

ggplot(dclimbefleaout, aes(x=tmeanLO, color=year)) +
    geom_histogram() +
    facet_wrap(.~lat, scales="free") +
    xlab(expression(paste("60 d mean temperature (", degree, "C)"), sep="")) +
    theme_minimal()


##
## Now try some simulations using the real climate data ...
fstar <- 200

dhere <- d[which(d$lat.lon %in% lookuppepid$lat.lon),]

simleafout <- data.frame(lat=numeric(),
                   lon=numeric(),
                   lat.lon=character(),
                   year=numeric(),
                   leafoutdoy=numeric(),
                   tmeanLO=numeric())

for(asite in unique(dhere$lat.lon)){
    onesitedf <- dhere[which(dhere$lat.lon==asite),]
    for(oneyear in unique(onesitedf$year)){
        oneyeardf <-  onesitedf[which(onesitedf$year==oneyear),]
        oneyeardf$tmeanforgdd <- ifelse(oneyeardf$tmean>0, oneyeardf$tmean, 0)
        leafoutdoy <-  min(which(cumsum(oneyeardf$tmeanforgdd) > fstar))
        climwindowsim <- oneyeardf[which(oneyeardf$doy<leafoutdoy & oneyeardf$doy>(leafoutdoy-windowlength)),]
        simloadd <- data.frame(
                   lat=oneyeardf$lat[1],
                   lon=oneyeardf$lon[1],
                   lat.lon=asite,
                   year=oneyear,
                   leafoutdoy=leafoutdoy,
                   tmeanLO=mean(climwindowsim$tmean, na.rm=TRUE))
        simleafout <- rbind(simleafout, simloadd)
    }
}

simleafoutwmean <- merge(dspcwindowmean, simleafout, by=c("lat", "lon", "lat.lon", "year"))
simleafoutwmean$tmeanplus10 <- simleafoutwmean$tmean+10
simleafoutwmean$tmeanLOplus10 <- simleafoutwmean$tmeanLO+10

simstatz <- data.frame(lat=numeric(),
                   lon=numeric(),
                   lat.lon=character(),
                   lmslopeSpcWindow=numeric(),
                   logslopeSpcWindow=numeric(),
                   lmslopeTLOd=numeric(),
                   logslopeTLOd=numeric())

for (site in unique(simleafoutwmean$lat.lon)){
    subby <- simleafoutwmean[which(simleafoutwmean$lat.lon==site),]
    statzadd <- data.frame(
                   lat=subby$lat[1],
                   lon=subby$lon[1],
                   lat.lon=site,
                   lmslopeSpcWindow=coef(lm(leafoutdoy~tmeanplus10, data=subby))[2],
                   logslopeSpcWindow=coef(lm(log(leafoutdoy)~log(tmeanplus10), data=subby))[2],
                   lmslopeTLOd =coef(lm(leafoutdoy~tmeanLOplus10, data=subby))[2],
                   logslopeTLOd=coef(lm(log(leafoutdoy)~log(tmeanLOplus10), data=subby))[2],
                   n=nrow(subby))
    simstatz <- rbind(simstatz, statzadd)
 }

pchhere <- 16
colz <- c("dodgerblue", "darkred", "darkorchid1")

pdf("pep_analyses/figures/loglinermodels_bylatitude_sims.pdf", width=8, height=8)
par(mfrow=c(2,2))
plot(lmslopeTLOd~lat, data=simstatz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean window before event",
     main="Linear -- Temp window before LO")
abline(lm(lmslopeTLOd~lat, data=simstatz), col=colz[1])

plot(logslopeTLOd~lat, data=simstatz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope w log: T mean window before event",
     main="Logged -- Temp window before LO")
abline(lm(logslopeTLOd~lat, data=simstatz), col=colz[2])

plot(lmslopeSpcWindow~lat, data=simstatz, col=colz[1], pch=pchhere,
     xlab="Latitude", ylab="Slope lm: T mean of Apr-June",
     main="Logged -- Temp window static")
abline(lm(lmslopeSpcWindow~lat, data=simstatz), col=colz[1])

plot(logslopeSpcWindow~lat, data=simstatz, col=colz[2], pch=pchhere,
     xlab="Latitude", ylab="Slope w log: T mean of Apr-June",
     main="Logged -- Temp window static")
abline(lm(logslopeSpcWindow~lat, data=simstatz), col=colz[2])
dev.off()


# Compare real and sim data ...
drealsim <- merge(dprepforstat, simleafoutwmean, by=c("lat", "lon", "lat.lon", "year"),
    all.x=TRUE, suffixes=c(".real", ".sim"))
statzrealsim <- merge(statz, simstatz, by=c("lat", "lon", "lat.lon"),
    all.x=TRUE, suffixes=c(".real", ".sim"))

par(mfrow=c(1,3))
plot(doy~leafoutdoy, data=drealsim) 
abline(0,1)
plot(tmean.real~tmean.sim, data=drealsim) # safety-check!
plot(tmeanLO.real~tmeanLO.sim, data=drealsim) 

par(mfrow=c(2,2))
plot(lmslopeSpcWindow.real~lmslopeSpcWindow.sim, data=statzrealsim)
abline(0,1)
plot(logslopeSpcWindow.real~logslopeSpcWindow.sim, data=statzrealsim) 
abline(0,1)
plot(lmslopeTLOd.real~lmslopeTLOd.sim, data=statzrealsim)
abline(0,1)
plot(logslopeTLOd.real~logslopeTLOd.sim, data=statzrealsim) 
abline(0,1)


plotleafoutrealsimsbylat <-  ggplot(drealsim, aes(x=leafoutdoy, y=doy, color=lat)) +
    geom_point() +
    xlab("Day of year: Simulated") +
    ylab(expression(paste("Day of year: real"), sep="")) +
    theme_minimal()

plotleafoutrealsimsbyyr <-  ggplot(drealsim, aes(x=leafoutdoy, y=doy, color=year)) +
    geom_point() +
    facet_wrap(.~lat) +
    xlab("Day of year: Simulated") +
    ylab(expression(paste("Day of year: real"), sep="")) +
    theme_minimal()


pdf("pep_analyses/figures/plotleafoutrealsimsbylat.pdf", width=8, height=8)
plotleafoutrealsimsbylat
dev.off()

pdf("pep_analyses/figures/plotleafoutrealsimsbyyr.pdf", width=8, height=8)
plotleafoutrealsimsbyyr
dev.off()



if(FALSE){
##
## Sims with sine curve
# Set up time periods and noise
dayz <- 366
yearz <- 50

sigma_y <- 5
x <- (1:(dayz*yearz))

# Set up one sine curve
amplitude =  20
center = 60
width =  dayz/2
offyosets <- c(0, 5, 10, 15, 20, 30)

# Set up for leafout
threshold <- 200
biozero <- -2 # above what temperature do plants accumulate

# The below does offset x year more efficiently
dfleafy <- data.frame(leafout_date = numeric(0),
                   cum_temp  = numeric(0),
                   mean_temp = numeric(0),
                   threemon_temp = numeric(0),
                   threshold = numeric(0),
                   offset = numeric(0))

for (offset in offyosets){
    for (i in 1:yearz){
        y = amplitude * sin(pi * (x - center) / width) + offset
        ybio <- y + rnorm(length(y), 0, sigma_y)
        ybio[which(ybio<biozero)] <- biozero
        leafout_date <- which.min(cumsum(ybio) < threshold)
        cum_temp <- sum(ybio[1:leafout_date])
        mean_temp <- mean(ybio[1:leafout_date])
        threemon_temp <- mean(ybio[30:120])
        dfleafy <- rbind(dfleafy, data.frame(leafout_date, cum_temp, mean_temp, threshold, offset))
        }
}

sapply(unique(dfleafy$offset), function(i) coef(lm(dfleafy$leafout_date[dfleafy$offset == i]~
    dfleafy$mean_temp[dfleafy$offset == i])))
}
