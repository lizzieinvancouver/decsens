## Started 25 January 2021 ##
## By Lizzie ##

## Some taken from analyseataglance.R 

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/projects/treegarden/decsens/analyses") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("hereliesboomboom")


library(ggplot2)

setwd("~/Documents/git/projects/vinmisc/vintages/analyses")

# climate data from Ben
temp <- read.csv("../../grapesdrought/WINELIZZIE/data/seas_temp_MJJ.onedeg.csv",
    header=TRUE)
colnames(temp)[1] <- "year"
prec <- read.csv("../../grapesdrought/WINELIZZIE/data/seas_prec_MJJ.onedeg.csv",
    header=TRUE)
colnames(prec)[1] <- "year"
pdsi <- read.csv("../../grapesdrought/WINELIZZIE/data/seas_pdsi_MJJ.onedeg.csv",
    header=TRUE)
colnames(pdsi)[1] <- "year"
daux <- read.csv("input/dauxdata.csv", header=TRUE, skip=2)
daux <- daux[,1:28]
names(daux)[names(daux)=="Abb."] <- "year"
ghd <- as.data.frame(daux)
ghd <- subset(ghd, select=c("year", "Bor", "Bur"))

envdata <- merge(ghd, temp, by="year", suffixes=c(".ghd", ".temp"))
envdata <- merge(envdata, prec, by="year", suffixes=c("", ".prec"))
envdata <- merge(envdata, pdsi, by="year", suffixes=c("", ".pdsi"))

## Setting up data to have same years 1980+ and before
subearly <- subset(envdata, year>1950 & year <1979)
sublate <- subset(envdata, year>1979)

nrow(subearly)
nrow(sublate)

summary(lm(Bur.ghd~Bur.temp, data=subearly))
summary(lm(Bur.ghd~Bur.temp, data=sublate))

summary(lm(log(Bur.ghd+60)~log(Bur.temp), data=subearly)) # adding 60 because the GHD can be negative
summary(lm(log(Bur.ghd+60)~log(Bur.temp), data=sublate))

summary(lm(Bor.ghd~Bor.temp, data=subearly))
summary(lm(Bor.ghd~Bor.temp, data=sublate))

summary(lm(log(Bor.ghd+60)~log(Bor.temp), data=subearly)) # adding 60 because the GHD can be negative
summary(lm(log(Bor.ghd+60)~log(Bor.temp), data=sublate))

# 20 years ...
subearly <- subset(envdata, year>1958 & year <1979)
sublate <- subset(envdata, year>1987)

## major axis
library(lmodel2)
# By Legendre https://cran.r-project.org/web/packages/lmodel2/lmodel2.pdf

MAearly <- lmodel2(Bur.ghd~Bur.temp, data=subearly, nperm=99)
MAlate <- lmodel2(Bur.ghd~Bur.temp, data=sublate, nperm=99)

wRMAearly <- lmodel2(Bur.ghd+60~Bur.temp, data=subearly, "relative", "relative", nperm=99)
wRMAlate <- lmodel2(Bur.ghd+60~Bur.temp, data=sublate, "relative", "relative", nperm=99)
