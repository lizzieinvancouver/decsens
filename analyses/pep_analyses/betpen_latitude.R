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

setwd("~/Documents/git/decsens/analyses/pep_analyses")
d<-read.csv("input/siteztotry.csv", header=TRUE)

rn<-brick("~/Desktop/Big Data Files/tn_0.25deg_reg_v16.0.nc", sep="")
rx<-brick("~/Desktop/Big Data Files/tx_0.25deg_reg_v16.0.nc", sep="")


period <- 1959:2011

d$lat.long <- paste(d$latnum, d$lonnum, sep="")
sites<-subset(d, select=c(latnum, lonnum, lat.long))
sites<-sites[!duplicated(sites$lat.long),]
sites$x<-sites$lonnum
sites$y<-sites$latnum
Coords<-subset(sites, select=c(x, y))
Coords <- na.omit(Coords)
nsites<-length(sites$lat.long)
tmin <- rn
tmax <- rx

points.min <- SpatialPoints(Coords, proj4string = rn@crs)
points.max <- SpatialPoints(Coords, proj4string = rx@crs)

yearsinclim<-as.numeric(format(as.Date(names(tmin),format="X%Y.%m.%d"),"%Y"))
yearsinperiod<-which(yearsinclim%in%period)
climsubmin<-subset(tmin,yearsinperiod)
climsubmax<-subset(tmax,yearsinperiod)

## subset climate days
monthsinclim<-as.numeric(format(as.Date(names(climsubmin),format="X%Y.%m.%d"),"%m"))
mstmonths<-c(3:4)
monthsinmst<-which(monthsinclim%in%mstmonths)
mstsubmin<-subset(climsubmin,monthsinmst)
mstsubmax<-subset(climsubmax,monthsinmst)

valuesmin <- raster::extract(mstsubmin,points.min)
valuesmax <- raster::extract(mstsubmax,points.max)

dclimmin <- cbind.data.frame(coordinates(points.min),valuesmin)
dclimmax <- cbind.data.frame(coordinates(points.max),valuesmax)

require(reshape2)
dxmin<-melt(dclimmin, id.vars=c("x","y"))
dxmax<-melt(dclimmax, id.vars=c("x","y"))

dxmin$lon <- dxmin$x
dxmin$lat <- dxmin$y
dxmin$date <- dxmin$variable
dxmin$Tmin <- dxmin$value

dxmax$lon <- dxmax$x
dxmax$lat <- dxmax$y
dxmax$date <- dxmax$variable
dxmax$Tmax <- dxmax$value

dx <- data.frame(lat=dxmin$lat, long=dxmin$lon, date=dxmin$date, tmin=dxmin$Tmin, tmax=dxmax$Tmax)
dx$Tavg <- (dx$tmin+dx$tmax)/2

dx$date<-substr(dx$date, 2,11)
dx$Date<- gsub("[.]", "-", dx$date)

dx$tmin <- NULL
dx$tmax <- NULL
dx$date<-NULL

dx$year<-as.numeric(substr(dx$Date, 0, 4))
dx$lat.long<-paste(dx$lat, dx$long)
dx$month <- substr(dx$Date, 6,7)
dx$doy <- as.numeric(strftime(dx$Date, format = "%j"))

### Now, let's vary pre-season length. We'll add 30, 45 and 60 days
dx$mat60<-ave(dx$Tavg, dx$year, dx$lat.long)

dailyandmst <- na.omit(dx)


write.csv(dailyandmst, file="~/Documents/git/decsens/analyses/pep_analyses/output/latitude_betpen_dailyandmst.csv", row.names = FALSE)

