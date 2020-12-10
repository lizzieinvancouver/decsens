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
setwd("~/Documents/git/decsens/analyses/pep_analyses")

###### FOR JUST ANALYZING RESULTS JUMP TO LINE 75 ########

# Get budburst data across 45 sites for BETPEN
# Betula puendula data from PEP (both have has GDD from 1 Jan to leafout)
# bp has mat from March 1st to June 1st and mat.lo is 30 days before leafout (uses tg -- aka mean -- data from E-OBS)
# bpalt is similar, but calculated uses txtm -- aka min and max (and we caculate the mean ourselves from those values) -- data from E-OBS) ... we don't use this currently 
d<-read.csv("input/pep_betpen_all.csv", header=TRUE)
#d<-read.csv("/n/wolkovich_lab/Lab/Cat/pep_betpen_all.csv", header=TRUE)

df<-d%>%
  filter(BBCH==11)%>%
  filter(YEAR>=1950 & YEAR<=2016)%>%
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
allpeps <- df[(df$year>=1951 & df$year<=2016),]

allpeps$cc<-NA
allpeps$cc<-ifelse(allpeps$year>1950 & allpeps$year<=1983, "1950-1983", allpeps$cc)
allpeps$cc<-ifelse(allpeps$year>1983 & allpeps$year<=2016, "1984-2016", allpeps$cc)
#allpeps$cc<-ifelse(allpeps$year>1970 & allpeps$year<=1990, "1970-1990", allpeps$cc)
allpeps$num.years<-ave(allpeps$year, allpeps$lat.long, FUN=length)
mostdata<-allpeps[(allpeps$num.years>=65),]
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

bbswpre <- bbsw[(bbsw$Year>1950 & bbsw$Year<=1983),]
bbswpost <- bbsw[(bbsw$Year>1983 & bbsw$Year<=2016),]
#bbswmid <- bbsw[(bbsw$Year>1970 & bbsw$Year<=1990),]

#bbswtest <- bbswpre[(bbswpre$spatial==1),]

### Now get the climate data for 45 sites for BETPEN (from betpen_climate_slidingwin.R)
#climatedatapre <- read.csv("/n/wolkovich_lab/Lab/Cat/bp_climatedatapre.csv")
climatedatapre <- read.csv("output/bp_climatedatapre.csv")
#climatedatapost <- read.csv("/n/wolkovich_lab/Lab/Cat/bp_climatedatapost.csv")
climatedatapost <- read.csv("output/bp_climatedatapost.csv")
#climatedatamid <- read.csv("/n/wolkovich_lab/Lab/Cat/bp_climatedatamid.csv")
#climatedatamid <- read.csv("output/bp_climatedatamid.csv")
  
source("simmonds_slidingwin/Run_SW.R")
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
write.csv(Results_SWRpre[[2]], file="simmonds_slidingwin/output/results_swapre_bp_mayref.csv")
write.csv(Results_SWRpre[[1]], file="simmonds_slidingwin/output/sumstats_swapre_bp_mayref.csv")

### Now checking Simmond's sliding window approach:
#refday <- c(01, 05) ### results in folders are from a ref day of 01-03, I think this new ref day is more appropriate for PEP leafout data - to rerun
datafile <- bbswpost
climate <- climatedatapost
climate$X <- NA ### needed in order to run... 

Results_SWRpost <- run_SW(absolute=TRUE, datafile, climate, refday) ## takes a long time to run
write.csv(Results_SWRpost[[2]], file="simmonds_slidingwin/output/results_swapost_bp_mayref.csv")
write.csv(Results_SWRpost[[1]], file="simmonds_slidingwin/output/sumstats_swapost_bp_mayref.csv")

if(FALSE){
### Now checking Simmond's sliding window approach:
#refday <- c(01, 04) ### results in folders are from a ref day of 01-03, I think this new ref day is more appropriate for PEP leafout data - to rerun
datafile <- bbswmid
climate <- climatedatamid
climate$X <- NA ### needed in order to run... 

Results_SWRmid <- run_SW(absolute=TRUE, datafile, climate, refday) ## takes a long time to run
write.csv(Results_SWRmid[[2]], file="/n/wolkovich_lab/Lab/Cat/results_swamid_bp_mayref.csv")
write.csv(Results_SWRmid[[1]], file="/n/wolkovich_lab/Lab/Cat/sumstats_swamid_bp_mayref.csv")
}

if(FALSE){
## Get data and parameters for prediction
source('pep_sims/simmonds_slidingwin/Params_SW.R')

# extract parameters for complete dataset
Parameters_SWRpost <- get_params_SW(Results_SWRpost, bbswpost$bb_mean, "complete", type = "Params")
Parameters_SWRpre <- read.csv("pep_sims/simmonds_slidingwin/output/parameters_swapre_mayref_bp.csv")
# SAVE
write.csv(Parameters_SWRpost, "pep_sims/simmonds_slidingwin/output/parameters_swapost_mayref_bp.csv", row.names=T)


######################################################################
#### Now we can analyze the results from the SWA: ref = May 31st #####
######################################################################

### Now let's check out the data
swapre <- read.csv("simmonds_slidingwin/output/results_swapre_bp_mayref.csv")
#swapre <- Results_SWRpre[[2]]
swapost <- read.csv("simmonds_slidingwin/output/results_swapost_bp_mayref.csv")
#swapost <- Results_SWRpost[[2]]

swapre_stat <- read.csv("simmonds_slidingwin/output/sumstats_swapre_bp_mayref.csv")
#swapre_stat <- Results_SWRpre[[1]]
swapost_stat <- read.csv("simmonds_slidingwin/output/sumstats_swapost_bp_mayref.csv")
#swapost_stat <- Results_SWRpost[[1]]

if(FALSE){
### Alright, now we have to convert to Kelvin again...
#swapre$tempk <- swapre$climate + 273.15
#swapost$tempk <- swapost$climate + 273.15

estprecc <- lm(yvar~climate, data=swapre)
estpostcc <- lm(yvar~climate, data=swapost)

diffbefore.after <- coef(estprecc)[2]-coef(estpostcc)[2]

estprecc.log <- lm(log(yvar)~log(climate), data=swapre)
estpostcc.log <- lm(log(yvar)~log(climate), data=swapost)

logdiffbefore.after <- coef(estprecc.log)[2]-coef(estpostcc.log)[2]

##### Now let's try and compare the model output from the SWA to mean, sd, and variance
mean(swapre$yvar) ### 113.8089
mean(swapost$yvar) ## 106.3356

mean(swapre$climate) ## 7.861508  
mean(swapost$climate) ## 5.675819 

sd(swapre$yvar) ## 11.58426
sd(swapost$yvar) ## 7.950618

sd(swapre$climate) ## 1.479552
sd(swapost$climate) ## 1.403269

var(swapre$yvar) ## 134.195
var(swapost$yvar) ## 63.21232

var(swapre$climate) ## 2.189074
var(swapost$climate) ## 1.969164

## Pre CC SWA results
#   X WindowOpen WindowClose Variable     Int      EST       SE        R2
#1 1        -22          52        1 157.0947 -5.50605 0.262991 0.4934151

## Post CC SWA results
#   WindowOpen WindowClose Variable   Int       EST        SE        R2
#1        -34          22        1 126.703 -3.588456 0.2071495 0.3998027
}

######################################################################
## Now using code from `bb_analysis/PEP_climate/comparetopepsims.R` ##
######################################################################
swapre$site <- rep(1:45)
swapost$site <- rep(1:45)

swa <- rbind(swapre, swapost)
swa$cc <- rep(c("pre", "post"), each=450)

swaest <- data.frame(site=numeric(), cc=character(), meanmat=numeric(), varmat=numeric(),  
                    sdmat=numeric(), 
                    matslope=numeric(), matslopese=numeric(), meanmatlo=numeric(), 
                    matslopelog=numeric(), matslopelogse=numeric())

sitez <- unique(swa$site)

for(i in c(1:length(sitez))){ # i <- 1
  subby <- subset(swa, site==sitez[i])
  for(ccstate in c(1:2)){
    subbycc <- subset(subby, cc==unique(swa$cc)[ccstate])
    meanmat <- mean(subbycc$climate, na.rm=TRUE)
    varmat <- var(subbycc$climate, na.rm=TRUE)
    sdmat <- sd(subbycc$climate, na.rm=TRUE)
    meanlo <- mean(subbycc$yvar, na.rm=TRUE)
    varlo <- var(subbycc$yvar, na.rm=TRUE)
    sdlo <- sd(subbycc$yvar, na.rm=TRUE)
    lmmat <- lm(yvar~climate, data=subbycc)
    lmmatse <- summary(lmmat)$coef[2,2]
    lmmatlog <- lm(log(yvar)~log(climate), data=subbycc)
    lmmatlogse <- summary(lmmatlog)$coef[2,2]
    swaestadd <- data.frame(siteslist=sitez[i], cc=unique(swa$cc)[ccstate], meanmat=meanmat, 
                           varmat=varmat, sdmat=sdmat, matslope=coef(lmmat)["climate"], matslopese=lmmatse, 
                           matslopelog=coef(lmmatlog)["log(climate)"], matslopelogse=lmmatlogse)
    swaest <- rbind(swaest, swaestadd)
  }
}    

swaest$matsloplog_exp <- exp(swaest$matslopelog)

write.csv(swaest, file=("output/swaestimates_withlog.csv"), row.names = FALSE)

#cc            meanbb     varbb      sdbb   meanclim   varclim   sdclim climslope   climslopese
#1 post       106.3356  45.47136  6.510241 14.49351  1.391118 1.110196 -2.048819  2.111988
#2  pre       113.8089 132.02667 11.087655 10.81471 10.969837 2.826528 -2.683118  1.468094


#   cc   meanbb    varbb  meanclim   varclim  climslope
#1 post 4.759305 27.06196 1.436775  1.062947 2.605993
#2  pre 3.926242 67.44772 1.637718 12.125168 1.923587



#### Compared to Lizzie's approach to BETPEN PEP data:
#          cc  meanmat   varmat    sdmat meanmatlo varmatlo  sdmatlo   meanlo     varlo     sdlo meanutah  meangdd  matslope matslopese
# 1950-1960 5.365163 3.005094 1.731358  6.814883 1.363054 1.086849 113.8089 110.51111 10.25803 2246.987 68.70881 -4.534630   1.258845
# 2000-2010 6.450939 1.251629 1.111780  6.615273 1.431603 1.152353 106.3356  46.95728  6.57374 2235.493 61.50754 -3.611025   1.579758



#### But when we compare to sliding window approach: with March 1st ref window
Parameters_SWRpre
#     WindowOpen WindowClose Variable      Int       EST        SE        R2
#1       -239        -238        1     136.3565 -2.084899 0.1192322 0.4043209

Parameters_SWRpost
#      WindowOpen WindowClose Variable      Int       EST        SE        R2
#1       -184        -176        2     143.0215 -2.531195 0.1695477 0.3307278
}


