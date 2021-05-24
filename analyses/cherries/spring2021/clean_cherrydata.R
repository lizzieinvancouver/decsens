### Clean up Cherry data 
# Started 24 May 2021 by Cat

## Housekeeping
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
library(tidyr)

setwd("~/Documents/git/decsens/analyses/cherries/spring2021/")


#### Start with the NOAA data...
noaa <- read.csv("input/swisscherry_noaa.csv")

noaa <- as.data.frame(noaa[-c(1:139),])
colnames(noaa) <- noaa[1,]

noaa_clean <- separate(noaa, `YEAR       OBS  OBSCORR    RECON     TT24`, c("year", "flo", "flo_corr", "flo_recon", "meantemp_febtoapr"))

noaa_clean <- noaa_clean[-1, ]
write.csv(noaa_clean, file="output/clean_cherry_noaa.csv", row.names=FALSE)


### Next with the PEP725 data for Germany and Norway...
ger <- read.csv("input/pruavi_germany_pep.csv")
gersites <- read.csv("input/PEP725_DE_stations.csv")
nor <- read.csv("input/pruavi_norway_pep.csv")
norsites <- read.csv("input/PEP725_NO_stations.csv")

ger.clean <- separate(ger, `PEP_ID.BBCH.YEAR.DAY`, c("pepid", "bbch", "year", "doy"))
ger.clean <- ger.clean[(ger.clean$bbch==60),]

gersites$PEP_ID.National_ID.LON.LAT.ALT.NAME <- rownames(gersites)
rownames(gersites) <- 1:nrow(gersites)

colnames(gersites) <- "PEP_ID;National_ID;LON;LAT;ALT;NAME"
gersites.clean <- separate(gersites, `PEP_ID;National_ID;LON;LAT;ALT;NAME`, c("pepid", "natid", "lon", "lat", "alt", "name"))
gersites.clean <- subset(gersites.clean, select=c("pepid", "lon", "lat"))
gersites.clean <-gersites.clean[!duplicated(gersites.clean),]

ger.clean <- full_join(ger.clean, gersites.clean)



