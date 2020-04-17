## Started 17 April 2020 ##
## By Lizzie, but not going well ##

## Following betpen_chillandgdd_tntx_forsims_odyssey.R some ... ##
## and my climate notes ... ##

## Get leafout date and daily climate (say, Jan 1 to end of Apr for now) for full year range of
# good BETPEN PEP sites ##

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
require(dplyr)
require(tidyr)
require(ggplot2)
require(lubridate)
require(raster)

setwd("~/Documents/git/projects/treegarden/decsens/analyses/pep_analyses/")
d <- read.csv("input/pep_betpen_all.csv")

df<-d%>%
  filter(BBCH==11)%>%
  filter(YEAR>=1950 & YEAR<=2010)%>%
  dplyr::select(YEAR, DAY, BBCH, PEP_ID, LAT, LON, species)%>%
  rename(year=YEAR)%>%
  rename(lo=DAY)%>%
  rename(lat=LAT)%>%
  rename(lon=LON)
df <- dplyr::select(df, year, PEP_ID, lat, lon, lo)

df <- df[!duplicated(df),]

x <- paste(df$year, df$lo)
df$date <-as.Date(strptime(x, format="%Y %j"))
df$Date <- as.character(df$date)
df$lat.lon <- paste(df$lat, df$lon)
allpeps <- df[(df$year>=1951 & df$year<=2010),]
allpeps$num.years<-ave(allpeps$year, allpeps$PEP_ID, FUN=length)
mostdata <- allpeps[(allpeps$num.years>=60),]

lositeyear <- subset(mostdata, select=c("lo", "lat", "lon", "lat.lon", "year"))
sitez <- subset(lositeyear, select=c("lon","lat"))
getsitez <- sitez[!duplicated(sitez$lat),] # bad, should add lon, but works for these sites

rg <- brick("/Volumes/OrangeFiend/climate/eobs/tg_ens_mean_0.25deg_reg_v20.0e.nc", sep="")

# extract data for relevant sites from Jan - May for now .. 
goober <- raster::extract(rg, getsitez) # this seems bad
