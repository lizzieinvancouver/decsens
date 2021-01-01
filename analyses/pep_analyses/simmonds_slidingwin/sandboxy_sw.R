## Started 30 December 2020 ##
## By Lizzie ##

## Playing around with sliding window stuffs ##
## This file runs  the sliding window from climwin for one site for Betpen (PEP) and ...
# makes some notes on what happens in Run_SW f(x) ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(climwin)
require(lubridate)
require(dplyr)
require(tidyr)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/treegarden/decsens/analyses/pep_analyses") 
} else setwd("~/Documents/git/decsens/analyses/pep_analyses")

## Below copied from bp_sw_simmonds.R ##
d <- read.csv("input/pep_betpen_all.csv", header=TRUE)

df <- d%>%
  filter(BBCH==11)%>%
  filter(YEAR>=1950 & YEAR<=2016)%>%
  dplyr::select(YEAR, DAY, BBCH, PEP_ID, LAT, LON, species)%>%
  rename(year=YEAR)%>%
  rename(lo=DAY)%>%
  rename(lat=LAT)%>%
  rename(long=LON)
## Hmm... can we sequence from budburst to leafout to find the number of freezes between?
df <- dplyr::select(df, year, PEP_ID, lat, long, lo)

df <- df[!duplicated(df),]

x <- paste(df$year, df$lo)
df$date <- as.Date(strptime(x, format="%Y %j"))
df$Date <-  as.character(df$date)
df$lat.long  <-  paste(df$lat, df$long)
allpeps  <-  df[(df$year>=1951 & df$year<=2016),]

allpeps$cc <- NA
allpeps$cc <- ifelse(allpeps$year>1950 & allpeps$year<=1983, "1950-1983", allpeps$cc)
allpeps$cc <- ifelse(allpeps$year>1983 & allpeps$year<=2016, "1984-2016", allpeps$cc)
#allpeps$cc <- ifelse(allpeps$year>1970 & allpeps$year<=1990, "1970-1990", allpeps$cc)
allpeps$num.years <- ave(allpeps$year, allpeps$lat.long, FUN=length)
mostdata <- allpeps[(allpeps$num.years>=65),]
tt <- as.data.frame(table(mostdata$cc, mostdata$lat.long))
tt <- tt[!(tt$Freq==0),]
bestsites <- as.data.frame(table(tt$Var2))
bestsites <- bestsites[(bestsites$Freq>1),]
bestsites  <-  bestsites$Var1

allpeps.subset <- mostdata[(mostdata$lat.long %in% bestsites),]

sites <- subset(allpeps.subset, select=c(lat, long, lat.long))
sites <- sites[!duplicated(sites$lat.long),]
badsites <- c("54.5 11.1", "49.7667 11.55", "47.8 11.0167") 
sites <- sites[!(sites$lat.long%in%badsites),]
sites$x <- sites$long
sites$y <- sites$lat
nsites <- length(sites$lat.long)
sites$siteslist <- 1:nsites

lositeyear <- subset(allpeps.subset, select=c("lo", "lat", "long", "lat.long", "year"))
lositeyear <- lositeyear[!duplicated(lositeyear),]
lositeyear <- left_join(lositeyear, sites)
lositeyear <- na.omit(lositeyear)
  
bbsw <- subset(lositeyear, select=c("year", "lo", "siteslist"))
bbsw$bb_date <- as.Date(bbsw$lo, origin=paste0(bbsw$year, "-01-01"))
bbsw$bb_date <- as.character(bbsw$bb_date)
bbsw$doy95 <- bbsw$lo - 4

bbsw <- subset(bbsw, select=c("year", "bb_date", "lo", "doy95", "siteslist"))
colnames(bbsw) <- c("Year", "bb_date", "bb_mean", "doy95", "spatial") # these are the specific columns it needs: year with a capital Y, bb_date (an actual date), bb_mean is leafout DOY (same as date but in doy), doy95 if leafout minus 4 (??), and spatial is different sites.
bbsw$bb_date <- as.character(bbsw$bb_date)

bbswpre <- bbsw[(bbsw$Year>1950 & bbsw$Year<=1983),]
bbswpost <- bbsw[(bbsw$Year>1983 & bbsw$Year<=2016),]

climatedatapre <- read.csv("output/bp_climatedatapre.csv")
climatedatapost <- read.csv("output/bp_climatedatapost.csv")


## Below adapted from Run_SW.R ##
datafile <- subset(bbswpre, spatial==1)
climate <- subset(climatedatapre, spatial==1) 
datafile <- as.data.frame(datafile)
datafile$bb_date <- strptime(datafile$bb_date, format = "%Y-%m-%d")
climate$date <- strptime(climate$date, format="%Y-%m-%d")

abspreBPsite1 <- slidingwin(
   baseline = lm(datafile$bb_mean ~ 1), # baseline model.
   xvar=list(climate$temp), 
   cdate = climate$date, 
   bdate = datafile$bb_date,
   cmissing=FALSE, cinterval="day", type="absolute",
   spatial=list(as.factor(datafile$spatial), as.factor(climate$spatial)),
    range=c(366,0), refday = c(7, 4), # 99 earliest day 
    stat=c("mean", "min", "max", "slope"),
    func="lin")

abspreBPsite1[[1]]$BestModel
abspreBPsite1[[2]]$BestModel

# You can do the below to then run some of Run_SW stuff ...
results <- abspreBPsite1

# Okay! So the results are given as:
# Best model is selected of the one's climwin says are best across 4 'aggregate statistics' (taken from Simmons et al. 2019, "The aggregate statistic used in this method is typically a sum, mean (e.g. of the daily mean, minimum, or maximum temperature), minimum (e.g. the minimum mean daily temperature reached in a focal window) or maximum (e.g. the max- imum mean daily temperature reached in a focal window) envi- ronmental value, or the slope of environmental change across the window (the gradient of a linear model of daily mean temperature against date within the focal window).") They are referenced this way in th package:
  # 1 = mean
  # 2 = max
  # 3 = min
  # 4 = slope
# WindowOpen (relative to refday)
# WindowClose (relative to refday)
# which statistic (1-4)
# model intercept
# model slope
# error of slope (I believe) 
# adj R squared

# So, for example, for abspreBPsite1
# highest R2 is slope (though it's 0.398 versus 0.373 for sum(max))
# window open 65 days before April 4 and closes 62 days before
