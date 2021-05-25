### Clean up Cherry data 
# Started 24 May 2021 by Cat

## Housekeeping
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
library(tidyr)
library(lubridate)

setwd("~/Documents/git/decsens/analyses/cherries/spring2021/")


#### Start with the NOAA data...
noaa <- read.csv("input/swisscherry_noaa.csv")

noaa <- as.data.frame(noaa[-c(1:139),])
colnames(noaa) <- noaa[1,]

noaa_clean <- separate(noaa, `YEAR       OBS  OBSCORR    RECON     TT24`, c("year", "flo", "flo_corr", "flo_recon", "meantemp_febtoapr"))

noaa_clean <- noaa_clean[-1, ]
write.csv(noaa_clean, file="output/clean_cherry_noaa.csv", row.names=FALSE)


#### Start with the NOAA data...
swiss <- read.csv("input/mateoswiss.csv")

colnames(swiss)

swiss_clean <- separate(swiss, `stn.time.mprua60d.mprua65d`, sep=";", c("station", "year", "flo_start", "flowering"))

swiss_clean$flo_start <- as.Date(swiss_clean$flo_start, "%Y%m%d")
swiss_clean$flowering <- as.Date(swiss_clean$flowering, "%Y%m%d")

swiss_clean$flo_start <- yday(swiss_clean$flo_start)
swiss_clean$flowering <- yday(swiss_clean$flowering)

write.csv(swiss_clean, file="output/clean_cherry_swiss.csv", row.names=FALSE)


### Finally with the PEP725 data for Germany and Norway...
ger <- read.csv("input/pruavi_germany_pep.csv")
gersites <- read.csv("input/PEP725_DE_stations.csv")

ger.clean <- separate(ger, `PEP_ID.BBCH.YEAR.DAY`, c("pepid", "bbch", "year", "doy"))
ger.clean <- ger.clean[(ger.clean$bbch==60),]

gersites$PEP_ID.National_ID.LON.LAT.ALT.NAME <- rownames(gersites)
rownames(gersites) <- 1:nrow(gersites)

colnames(gersites) <- "PEP_ID;National_ID;LON;LAT;ALT;NAME"
gersites.clean <- separate(gersites, `PEP_ID;National_ID;LON;LAT;ALT;NAME`, sep=";", c("pepid", "natid", "lon", "lat", "alt", "name"))
gersites.clean <- subset(gersites.clean, select=c("pepid", "lon", "lat"))
gersites.clean <-gersites.clean[!duplicated(gersites.clean),]

ger.clean <- left_join(ger.clean, gersites.clean)


nor <- read.csv("input/pruavi_norway_pep.csv")
norsites <- read.csv("input/PEP725_NO_stations.csv")


nor.clean <- separate(nor, `PEP_ID.BBCH.YEAR.DAY`, c("pepid", "bbch", "year", "doy"))
nor.clean <- nor.clean[(nor.clean$bbch==60),]

colnames(norsites) <- "PEP_ID;National_ID;LON;LAT;ALT;NAME"
norsites.clean <- separate(norsites, `PEP_ID;National_ID;LON;LAT;ALT;NAME`, sep=";", c("pepid", "natid", "lon", "lat", "alt", "name"))
norsites.clean <- subset(norsites.clean, select=c("pepid", "lon", "lat"))
norsites.clean <-norsites.clean[!duplicated(norsites.clean),]

nor.clean <- left_join(nor.clean, norsites.clean)


pepclean <- full_join(ger.clean, nor.clean)


write.csv(pepclean, file="output/clean_cherry_pep.csv", row.names=FALSE)
