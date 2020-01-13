## Started 13 Jan 2020
## By Cat - based off code by Lizzie in ospree repo ospree/analyses/bb_analysis/pep_sims/comparetopepsims.R

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

#Libraries
library(dplyr)
library(tidyr)

# Setting working directory. Add in your own path in an if statement for your file structure
setwd("~/Documents/git/decsens/analyses/pep_analyses") 

# get some data
# Betula puendula data from PEP (both have has GDD from 1 Jan to leafout)
# bp has mat from March 1st to June 1st and mat.lo is 30 days before leafout (uses tg -- aka mean -- data from E-OBS)
# bpalt is similar, but calculated uses txtm -- aka min and max (and we caculate the mean ourselves from those values) -- data from E-OBS) ... we don't use this currently 
bp <- read.csv("output/betpen_allchillsandgdds_45sites_mat_tntx_forsims.csv", header=TRUE)
bppre <- read.csv("output/bp_climatedatapre.csv")
bppost <- read.csv("output/bp_climatedatapost.csv")

# loop to extract some model estimates
# this takes mean for each time period then allows comparison acrosgs the two resulting values
bpest <- data.frame(siteslist=numeric(), cc=character(), meanmat=numeric(), varmat=numeric(),  
                    sdmat=numeric(), meanlo=numeric(), varlo=numeric(), sdlo=numeric(), meanutah=numeric(), meangdd=numeric(), 
                    matslope=numeric(), matslopese=numeric(), meanmatlo=numeric(), varmatlo=numeric(), sdmatlo=numeric())

sitez <- unique(bp$siteslist)

for(i in c(1:length(sitez))){ # i <- 1
  subby <- subset(bp, siteslist==sitez[i])
  for(ccstate in c(1:2)){
    subbycc <- subset(subby, cc==unique(bp$cc)[ccstate])
    meanmat <- mean(subbycc$mat, na.rm=TRUE)
    varmat <- var(subbycc$mat, na.rm=TRUE)
    sdmat <- sd(subbycc$mat, na.rm=TRUE)
    meanmatlo <- mean(subbycc$mat.lo, na.rm=TRUE)
    varmatlo <- var(subbycc$mat.lo, na.rm=TRUE)
    sdmatlo <- sd(subbycc$mat.lo, na.rm=TRUE)
    meanlo <- mean(subbycc$lo, na.rm=TRUE)
    varlo <- var(subbycc$lo, na.rm=TRUE)
    sdlo <- sd(subbycc$lo, na.rm=TRUE)
    meanutah <- mean(subbycc$chillutah, na.rm=TRUE)
    meangdd <- mean(subbycc$gdd, na.rm=TRUE)
    lmmat <- lm(lo~mat, data=subbycc)
    lmmatse <- summary(lmmat)$coef[2,2]
    bpestadd <- data.frame(siteslist=sitez[i], cc=unique(bp$cc)[ccstate], meanmat=meanmat, 
                           varmat=varmat, sdmat=sdmat, meanlo=meanlo, varlo=varlo, sdlo=sdlo, meanutah=meanutah, 
                           meangdd=meangdd, matslope=coef(lmmat)["mat"], matslopese=lmmatse, meanmatlo=meanmatlo,
                           varmatlo=varmatlo, sdmatlo=sdmatlo)
    bpest <- rbind(bpest, bpestadd)
  }
}    

meanhere <- aggregate(bpest[c("meanmat", "varmat", "sdmat", "meanmatlo", "varmatlo", "sdmatlo", "meanlo", "varlo", "sdlo", "meanutah", "meangdd",
                              "matslope", "matslopese")], bpest["cc"], FUN=mean)
sdhere <- aggregate(bpest[c("meanmat", "varmat", "meanmatlo", "varmatlo", "meanlo", "varlo", "meanutah", "meangdd", "matslope")],
                    bpest["cc"], FUN=sd)


#          cc  meanmat   varmat    sdmat meanmatlo varmatlo  sdmatlo   meanlo     varlo     sdlo meanutah  meangdd  matslope matslopese
# 1950-1960 5.365163 3.005094 1.731358  6.814883 1.363054 1.086849 113.8089 110.51111 10.25803 2246.987 68.70881 -4.534630   1.258845
# 2000-2010 6.450939 1.251629 1.111780  6.615273 1.431603 1.152353 106.3356  46.95728  6.57374 2235.493 61.50754 -3.611025   1.579758


## Also get the difference for each site across two time periods
# This is to compare to sims better

bpest.sitediffs <- data.frame(siteslist=numeric(), matdiff=numeric(), matlodiff=numeric(), diffslope=numeric(),
                              varlodiff=numeric(), varlodiffper=numeric(), varmatdiffper=numeric())

for(i in c(1:length(sitez))){ # i <- 1
  subby <- subset(bpest, siteslist==sitez[i])
  precc <- subset(subby, cc=="1950-1960")
  postcc <- subset(subby, cc=="2000-2010")
  matdiff <- precc$meanmat-postcc$meanmat
  matlodiff <- precc$meanmatlo-postcc$meanmatlo
  diffslope <- precc$matslope-postcc$matslope
  varlodiff <- precc$varlo-postcc$varlo
  varlodiffper <- postcc$varlo/precc$varlo
  varmatdiffper <- postcc$varmat/precc$varmat
  bpest.sitediffs.add <- data.frame(siteslist=sitez[i], matdiff=matdiff,matlodiff=matlodiff, diffslope=diffslope,
                                    varlodiff=varlodiff, varlodiffper=varlodiffper, varmatdiffper=varmatdiffper)
  bpest.sitediffs <- rbind(bpest.sitediffs, bpest.sitediffs.add)
}

bpest.sitediffs$daysperC <- bpest.sitediffs$diffslope/bpest.sitediffs$matdiff


bpclimall <- full_join(bppre, bppost)
bpclimall <- subset(bpclimall, select=c(yday, year, temp, spatial))
bpclimall$cc <- ifelse(bpclimall$year>=1950 & bpclimall$year<=1960, "1950-1960", "2000-2010")
names(bpclimall) <- c("doy", "year", "temp", "siteslist", "cc")
bpclimall <- bpclimall[!duplicated(bpclimall),]
nonyrs <- c(1950, 1961, 2000, 2011)
bpclimall <- bpclimall[!(bpclimall$year%in%nonyrs),]
#bpclimall <- bpclimall[!(bpclimall$doy<=60),]

lodate <- subset(bp, select=c("year", "siteslist", "cc", "lo"))
lodate <- lodate[!duplicated(lodate),]

lopersite <- data.frame()
dailytemps <- data.frame()
for(i in length(sitez)){ #i=1
  for(j in 1951:2010){ #j=1951
    lopersite <- lodate[(lodate$siteslist==i & lodate$year==j),]
    lo <- as.numeric(lopersite$lo)
    dailytemps <- bpclimall[(bpclimall$doy<=lo),]
    
  }
}

bptempandbb <- full_join(dailytemps, lodate)
bptempandbb <- bptempandbb[!duplicated(bptempandbb),]
bptempandbb$tempk <- bptempandbb$temp + 273.15

#### Now calculate temperature sensitivites but with log
estprecc.log <- lm(log(lo)~log(tempk), data=subset(bptempandbb, cc=="1950-1960")) 
estpostcc.log <- lm(log(lo)~log(tempk), data=subset(bptempandbb, cc=="2000-2010"))

logdiffbefore.after <- coef(estprecc.log)[2]-coef(estpostcc.log)[2]
# negative means a decline in sensitivity AFTER climate change

#### Now calculate temperature sensitivites without log to compare
estprecc <- lm(lo~tempk, data=subset(bptempandbb, cc=="1950-1960")) 
estpostcc <- lm(lo~tempk, data=subset(bptempandbb, cc=="2000-2010"))

diffbefore.after <- coef(estprecc)[2]-coef(estpostcc)[2]
# negative means a decline in sensitivity AFTER climate change


########################
## Stuff in the supp ##
########################

# estimate change in variance
# take the means, then calculates the percent change between pre and post cc 
1-meanhere$varlo[2]/meanhere$varlo[1]
1-mean.pepsims$var.lo.postcc[1]/mean.pepsims$var.lo.precc[1]
1-median.pepsims$var.lo.postcc[1]/median.pepsims$var.lo.precc[1]

# in supp 'we estimated a decline in sensitivity'
mean(bpest.sitediffs$daysperC)
sd(bpest.sitediffs$daysperC)/sqrt(45) # days per C mean SE bp

# in supp 'given X warming'
mean(bpest.sitediffs$matdiff)
mean(bpest.sitediffs$matdiff)/sqrt(45)

# compare to pep sims (in figure, not text in supp currently)
mean.pepsims$diffbefore.after[1]
sd.pepsims$diffbefore.after[1]/(sqrt(45)) # days per C mean SE pepsims

# Utah units
meanhere$meanutah
sdhere$meanutah/sqrt(45)

# gdd and matlo
meanhere$meangdd
sdhere$meangdd/sqrt(45)
meanhere$meanmatlo
sdhere$meanmatlo/sqrt(45)

##############################
## End of stuff in the supp ##
##############################

# variance! here's what you get if you calculate the % change site-by-site then average
1 - mean(bpest.sitediffs$varlodiffper)
sd(bpest.sitediffs$varlodiffper)/sqrt(45)
1 - mean.pepsims$varlodiffper[1]
sd.pepsims$varlodiffper[1]/(sqrt(45))

# closer look at variance differences due to when we take the mean ...
# when we take the diff at each site, we exacerbate outliers, this explains (I think) the difference between the two ways we calculate variance
pepsims.1d <- subset(pepsims, degwarm==1)
hist(pepsims.1d$varlodiff)
hist(pepsims.1d$var.lo.precc)
hist(pepsims.1d$var.lo.postcc)

# other stuff ... (not using currently)
mean(bpest.sitediffs$diffslope)
sd(bpest.sitediffs$diffslope)/sqrt(45)
mean(bpest.sitediffs$matlodiff)
sd(bpest.sitediffs$matdiff)
sd(bpest.sitediffs$matlodiff)
sd(bpest.sitediffs$daysperC)

##############
## Plotting ##
##############

cexhere <- 0.95
pdf(file.path("figures/peprealandsims.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(0.5,4.5), ylim=c(-3.1, -0.1),
     ylab=expression(paste("Change in estimated sensitivity (days/", degree, "C)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:4){
  pos.x <- mean.pepsims$degwarm[i]
  pos.y <- mean.pepsims$diffbefore.after[i]
  sehere <- sd.pepsims$diffbefore.after[i]/(sqrt(45))
  lines(x=rep(pos.x, 2), y=c(pos.y-sehere, pos.y+sehere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
}
realdat.diff <- mean(bpest.sitediffs$diffslope) 
points(abs(mean(bpest.sitediffs$matdiff)), realdat.diff, cex=cexhere, pch=17, col="salmon")
realdatse <- sd(bpest.sitediffs$diffslope)/sqrt(45)
lines(x=rep(abs(mean(bpest.sitediffs$matdiff)), 2), y=c(realdat.diff-realdatse, realdat.diff+realdatse),
      col="salmon")
# par(xpd=TRUE) # so I can plot legend outside
legend("topright", pch=c(17, 19), col=c("salmon", "darkblue"), legend=c("European data (PEP725)", "Simulations with constant cues"),
       cex=1, bty="n")
dev.off()
