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

burearlylm <- lm(Bur.ghd~Bur.temp, data=subearly)
burlatelm <- lm(Bur.ghd~Bur.temp, data=sublate)

burearlylog <- lm(log(Bur.ghd+60)~log(Bur.temp), data=subearly) # adding 60 because the GHD can be negative
burlatelog <-  lm(log(Bur.ghd+60)~log(Bur.temp), data=sublate)

bordearlylm <- lm(Bor.ghd~Bor.temp, data=subearly)
bordlatelm <- lm(Bor.ghd~Bor.temp, data=sublate)

bordearlylog <- lm(log(Bor.ghd+60)~log(Bor.temp), data=subearly) # adding 60 because the GHD can be negative
bordlatelog <- lm(log(Bor.ghd+60)~log(Bor.temp), data=sublate)

burearly <- data.frame(where="Burgundy",
                    when="Before",
                    meantemp=mean(subearly$Bur.temp, na.rm=TRUE),
                    lmslope=coef(burearlylm)[2],
                    lmslopeconfint11=confint(burearlylm,level=0.89)[2,1],
                    lmslopeconfint89=confint(burearlylm,level=0.89)[2,2],
                    logslope=coef(burearlylog)[2],
                    logslopeconfint11=confint(burearlylog,level=0.89)[2,1],
                    logslopeconfint89=confint(burearlylog,level=0.89)[2,2])

burlate <- data.frame(where="Burgundy",
                    when="After",
                    meantemp=mean(sublate$Bur.temp, na.rm=TRUE),
                    lmslope=coef(burlatelm)[2],
                    lmslopeconfint11=confint(burlatelm,level=0.89)[2,1],
                    lmslopeconfint89=confint(burlatelm,level=0.89)[2,2],
                    logslope=coef(burlatelog)[2],
                    logslopeconfint11=confint(burlatelog,level=0.89)[2,1],
                    logslopeconfint89=confint(burlatelog,level=0.89)[2,2])

bordearly <- data.frame(where="Bordeaux",
                    when="Before",
                    meantemp=mean(subearly$Bor.temp, na.rm=TRUE),
                    lmslope=coef(bordearlylm)[2],
                    lmslopeconfint11=confint(bordearlylm,level=0.89)[2,1],
                    lmslopeconfint89=confint(bordearlylm,level=0.89)[2,2],
                    logslope=coef(bordearlylog)[2],
                    logslopeconfint11=confint(bordearlylog,level=0.89)[2,1],
                    logslopeconfint89=confint(bordearlylog,level=0.89)[2,2])
bordlate <- data.frame(where="Bordeaux",
                    when="After",
                    meantemp=mean(sublate$Bor.temp, na.rm=TRUE),
                    lmslope=coef(bordlatelm)[2],
                    lmslopeconfint11=confint(bordlatelm,level=0.89)[2,1],
                    lmslopeconfint89=confint(bordlatelm,level=0.89)[2,2],
                    logslope=coef(bordlatelog)[2],
                    logslopeconfint11=confint(bordlatelog,level=0.89)[2,1],
                    logslopeconfint89=confint(bordlatelog,level=0.89)[2,2])

plotdat <- rbind(burearly, burlate, bordearly, bordlate)

if(FALSE){
cexhere <- 0.95
cextext <- 0.75
pdf(file.path("figures/winedata.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
for(i in 1:length(unique(plotdat$when))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slope[i]
  ciherelo <- mean.betpen.20yr$mat60slopeconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen.20yr$cc)[i], cex=cextext, col="darkblue")
  }
for(i in 1:length(unique(mean.betpen.20yr$cc))){
  pos.x <- tempdiffplot.20yr[i]
  pos.y <- mean.betpen.20yr$mat60slopelog[i]
  ciherelo <- mean.betpen.20yr$mat60slopelogconfint11[i]
  cihereup <- mean.betpen.20yr$mat60slopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=1, bty="n")
dev.off()
}

## Side note ... 
## major axis regression (Type II)
if(FALSE){
library(lmodel2)
# By Legendre https://cran.r-project.org/web/packages/lmodel2/lmodel2.pdf
MAearly <- lmodel2(Bur.ghd~Bur.temp, data=subearly, nperm=99)
MAlate <- lmodel2(Bur.ghd~Bur.temp, data=sublate, nperm=99)

wRMAearly <- lmodel2(Bur.ghd+60~Bur.temp, data=subearly, "relative", "relative", nperm=99)
wRMAlate <- lmodel2(Bur.ghd+60~Bur.temp, data=sublate, "relative", "relative", nperm=99)
}
