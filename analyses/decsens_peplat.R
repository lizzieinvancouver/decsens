## Started 10 February 2021 ##
## By Lizzie, with help from Cat ##

## See also projects/misc/pep725/pep725check_decsens.R ##

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/projects/treegarden/decsens/analyses") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("hereliesboomboom")

d <- read.csv("pep_analyses/output/latitude_betpen_dailytemps.csv")
dm <- read.csv("pep_analyses/output/latitude_betpen_dailyandmst.csv")

d$date <- as.Date(d$Date, format=("%Y-%m-%d"))
d$doy <- format(d$date, "%j")

dsummdoy <- aggregate(d["tmean"], d[c("lat", "lon","lat.lon", "doy")], FUN=mean)

#
library(ggplot2)

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
