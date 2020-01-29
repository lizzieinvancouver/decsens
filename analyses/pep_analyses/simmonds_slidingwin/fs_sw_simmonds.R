# Started 29 October 2019 by Cat
## Using code from Emily Simmond's sliding window approach: https://github.com/emilygsimmonds/Cue_Identification

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(climwin)
require(lubridate)
require(dplyr)
require(tidyr)

# Setting working directory. Add in your own path in an if statement for your file structure
#setwd("~/Documents/git/decsens/analyses/pep_analyses")

###### FOR JUST ANALYZING RESULTS JUMP TO LINE 75 ########

# Get budburst data across 45 sites for BETPEN
# Betula puendula data from PEP (both have has GDD from 1 Jan to leafout)
# bp has mat from March 1st to June 1st and mat.lo is 30 days before leafout (uses tg -- aka mean -- data from E-OBS)
# bpalt is similar, but calculated uses txtm -- aka min and max (and we caculate the mean ourselves from those values) -- data from E-OBS) ... we don't use this currently 
d<-read.csv("/n/wolkovich_lab/Lab/Cat/pep_fagsyl_all.csv", header=TRUE)

df<-d%>%
  filter(BBCH==11)%>%
  filter(YEAR>=1950 & YEAR<=2010)%>%
  dplyr::select(YEAR, DAY, BBCH, PEP_ID, LAT, LON, species)%>%
  rename(year=YEAR)%>%
  rename(lo=DAY)%>%
  rename(lat=LAT)%>%
  rename(long=LON)
## Hmm... can we sequence from budburst to leafout to find the number of freezes between?
df<-dplyr::select(df, year, PEP_ID, lat, long, lo)

df<-df[!duplicated(df),]

x<-paste(df$year, df$lo)
df$date<-as.Date(strptime(x, format="%Y %j"))
df$Date<- as.character(df$date)
df$lat.long <- paste(df$lat, df$long)
allpeps <- df[(df$year>=1951 & df$year<=2010),]

allpeps$cc<-NA
allpeps$cc<-ifelse(allpeps$year>1950 & allpeps$year<=1970, "1950-1970", allpeps$cc)
allpeps$cc<-ifelse(allpeps$year>1990 & allpeps$year<=2010, "1990-2010", allpeps$cc)
allpeps$cc<-ifelse(allpeps$year>1970 & allpeps$year<=1990, "1970-1990", allpeps$cc)
allpeps$num.years<-ave(allpeps$year, allpeps$PEP_ID, FUN=length)
mostdata<-allpeps[(allpeps$num.years>=60),]
tt<-as.data.frame(table(mostdata$cc, mostdata$lat.long))
tt<-tt[!(tt$Freq==0),]
bestsites<-as.data.frame(table(tt$Var2))
bestsites<-bestsites[(bestsites$Freq>1),]
bestsites <- bestsites$Var1

allpeps.subset<-mostdata[(mostdata$lat.long %in% bestsites),]

sites<-subset(allpeps.subset, select=c(lat, long, lat.long))
sites<-sites[!duplicated(sites$lat.long),]
badsites<-c("54.5 11.1", "49.7667 11.55", "47.8 11.0167") 
sites<-sites[!(sites$lat.long%in%badsites),]
sites$x<-sites$long
sites$y<-sites$lat
nsites<-length(sites$lat.long)
sites$siteslist<-1:nsites

lositeyear <- subset(allpeps.subset, select=c("lo", "lat", "long", "lat.long", "year"))
lositeyear <- lositeyear[!duplicated(lositeyear),]
lositeyear <- left_join(lositeyear, sites)
lositeyear<-na.omit(lositeyear)

bbsw <- subset(lositeyear, select=c("year", "lo", "siteslist"))
bbsw$bb_date <- as.Date(bbsw$lo, origin=paste0(bbsw$year, "-01-01"))
bbsw$bb_date <- as.character(bbsw$bb_date)
bbsw$doy95 <- bbsw$lo - 4

bbsw <- subset(bbsw, select=c("year", "bb_date", "lo", "doy95", "siteslist"))
colnames(bbsw) <- c("Year", "bb_date", "bb_mean", "doy95", "spatial")
bbsw$bb_date <- as.character(bbsw$bb_date)

bbswpre <- bbsw[(bbsw$Year>1950 & bbsw$Year<=1970),]
bbswpost <- bbsw[(bbsw$Year>1990 & bbsw$Year<=2010),]
bbswmid <- bbsw[(bbsw$Year>1970 & bbsw$Year<=1990),]

#bbswtest <- bbswpre[(bbswpre$spatial==1),]

### Now get the climate data for 45 sites for BETPEN (from betpen_climate_slidingwin.R)
climatedatapre <- read.csv("/n/wolkovich_lab/Lab/Cat/fs_climatedatapre.csv")
climatedatapost <- read.csv("/n/wolkovich_lab/Lab/Cat/fs_climatedatapost.csv")
climatedatamid <- read.csv("/n/wolkovich_lab/Lab/Cat/fs_climatedatamid.csv")

source("/n/wolkovich_lab/Lab/Cat/Run_SW.R")
# refday = c(day, mon)
# climate is a datafile that must include col = temp
# datafile = biological data
# default = absolute but can also run relative
#run_SW <- function(absolute = TRUE, datafile, climate, refday)

### Now checking Simmond's sliding window approach:
refday <- c(01, 05) ### results in folders are from a ref day of 01-03, I think this new ref day is more appropriate for PEP leafout data - to rerun
datafile <- bbswpre
climate <- climatedatapre
climate$X <- NA ### needed in order to run... 

Results_SWRpre <- run_SW(absolute=TRUE, datafile, climate, refday) ## takes a long time to run
write.csv(Results_SWRpre[[2]], file="/n/wolkovich_lab/Lab/Cat/results_swapre_fs_mayref.csv")
write.csv(Results_SWRpre[[1]], file="/n/wolkovich_lab/Lab/Cat/sumstats_swapre_fs_mayref.csv")

### Now checking Simmond's sliding window approach:
refday <- c(01, 05) ### results in folders are from a ref day of 01-03, I think this new ref day is more appropriate for PEP leafout data - to rerun
datafile <- bbswpost
climate <- climatedatapost
climate$X <- NA ### needed in order to run... 

Results_SWRpost <- run_SW(absolute=TRUE, datafile, climate, refday) ## takes a long time to run
write.csv(Results_SWRpost[[2]], file="/n/wolkovich_lab/Lab/Cat/results_swapost_fs_mayref.csv")
write.csv(Results_SWRpost[[1]], file="/n/wolkovich_lab/Lab/Cat/sumstats_swapost_fs_mayref.csv")

### Now checking Simmond's sliding window approach:
refday <- c(01, 05) ### results in folders are from a ref day of 01-03, I think this new ref day is more appropriate for PEP leafout data - to rerun
datafile <- bbswmid
climate <- climatedatamid
climate$X <- NA ### needed in order to run... 

Results_SWRmid <- run_SW(absolute=TRUE, datafile, climate, refday) ## takes a long time to run
write.csv(Results_SWRmid[[2]], file="/n/wolkovich_lab/Lab/Cat/results_swamid_fs_mayref.csv")
write.csv(Results_SWRmid[[1]], file="/n/wolkovich_lab/Lab/Cat/sumstats_swamid_fs_mayref.csv")

if(FALSE){
## Get data and parameters for prediction
source('pep_sims/simmonds_slidingwin/Params_SW.R')

# extract parameters for complete dataset
Parameters_SWRpost <- get_params_SW(Results_SWRpost, bbswpost$bb_mean, "complete", type = "Params")
Parameters_SWRpre <- read.csv("pep_sims/simmonds_slidingwin/output/parameters_swapre_mayref_fs.csv")
# SAVE
write.csv(Parameters_SWRpost, "pep_sims/simmonds_slidingwin/output/parameters_swapost_mayref_fs.csv", row.names=T)
}