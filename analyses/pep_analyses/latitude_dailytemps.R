### Latitudinal Gradient test
### 8 Feb 2021 - Cat
## based off of betpen_chillandgdd_tntx_forsims.R but using Tmin and Tmax now to find Tmean

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
require(lubridate)
require(chillR)
require(raster)
library(ncdf4)

setwd("~/Documents/git/decsens/analyses/pep_analyses")
d<-read.csv("input/siteztotry.csv", header=TRUE)

d$lat.lon <- paste(d$latnum, d$lonnum, sep="")
getsitez <- d[!duplicated(d$lat.lon),] # bad, should add lon, but works for these sites - changed column to lat.lon and then address below
getsitez$x <- getsitez$lonnum
getsitez$y  <- getsitez$latnum
Coords <- subset(getsitez, select=c(x, y))
Coords <- na.omit(Coords)
nsites <- length(getsitez$lat.lon)

rg <- brick("~/Desktop/Big Data Files/tg_0.25deg_reg_v19.0.nc")

period <- c(1959:2011)

points <- SpatialPoints(Coords, proj4string = rg@crs)

yearsinclim <- as.numeric(format(as.Date(names(rg),format="X%Y.%m.%d"),"%Y"))
yearsinperiod <- which(yearsinclim%in%period)
climsub <- subset(rg,yearsinperiod) # takes a couple of minutes

## subset climate days
monthsinclim <- as.numeric(format(as.Date(names(climsub),format="X%Y.%m.%d"),"%m"))
dailytempmonths <- c(1:6)
monthsindailytemps <- which(monthsinclim%in%dailytempmonths)
dailytempssub <- subset(climsub,monthsindailytemps)


values <- raster::extract(dailytempssub,points) ## takes many minutes!

dclim <- cbind.data.frame(coordinates(points),values)

require(reshape2)
extclimdata <- melt(dclim, id.vars=c("x","y"))

extclimdata$lon <- extclimdata$x
extclimdata$lat <- extclimdata$y
extclimdata$date <- extclimdata$variable
extclimdata$Tmean <- extclimdata$value

dailyclimdata <- data.frame(lat=extclimdata$lat, lon=extclimdata$lon, date=extclimdata$date, tmean=extclimdata$Tmean)

dailyclimdata$date <- substr(dailyclimdata$date, 2,11)
dailyclimdata$Date <- gsub("[.]", "-", dailyclimdata$date)

dailyclimdata$date <- NULL

dailyclimdata$year <- as.numeric(substr(dailyclimdata$Date, 0, 4))
dailyclimdata$lat.lon <- paste(dailyclimdata$lat, dailyclimdata$lon)

dailyclimdata <- dailyclimdata[!duplicated(dailyclimdata),]

write.csv(dailyclimdata, file="~/Documents/git/decsens/analyses/pep_analyses/output/latitude_betpen_dailytemps.csv", row.names = FALSE)

