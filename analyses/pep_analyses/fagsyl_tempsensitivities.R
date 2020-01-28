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
# Fagus sylvatica data from PEP (both have has GDD from 1 Jan to leafout)
# fs has mat from March 1st to June 1st and mat.lo is 30 days before leafout (uses tg -- aka mean -- data from E-OBS)
# fsalt is similar, but calculated uses txtm -- aka min and max (and we caculate the mean ourselves from those values) -- data from E-OBS) ... we don't use this currently 
#fs <- read.csv("output/betpen_allchillsandgdds_45sites_mat_tntx_forsims.csv", header=TRUE)
#fspre <- read.csv("output/fs_climatedatapre.csv")
#fspost <- read.csv("output/fs_climatedatapost.csv")

fs <- read.csv("output/fagsyl_decsens_1950-1970_1990-2000.csv")

# loop to extract some model estimates
# this takes mean for each time period then allows comparison acrosgs the two resulting values
fsest <- data.frame(siteslist=numeric(), cc=character(), meanmat=numeric(), varmat=numeric(),  
                    sdmat=numeric(), meanlo=numeric(), varlo=numeric(), sdlo=numeric(), meanutah=numeric(), meangdd=numeric(), 
                    matslope=numeric(), matslopese=numeric(), matslopeconfint11=numeric(), matslopeconfint89=numeric(),
                    meanmatlo=numeric(), 
                    matslopelog=numeric(), matslopelogse=numeric(), matslopelogconfint11=numeric(), matslopelogconfint89=numeric(),
                    varmatlo=numeric(), sdmatlo=numeric())

sitez <- unique(fs$siteslist)

for(i in c(1:length(sitez))){ # i <- 1
  subby <- subset(fs, siteslist==sitez[i])
  for(ccstate in c(1:2)){ ## ccstate=1
    subbycc <- subset(subby, cc==unique(fs$cc)[ccstate])
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
    lmmatconfint11 <- confint(lmmat,level=89)[2,1]
    lmmatconfint89 <- confint(lmmat,level=89)[2,2]
    lmmatlog <- lm(log(lo)~log(mat), data=subbycc)
    lmmatlogse <- summary(lmmatlog)$coef[2,2]
    lmmatconfintlog11 <- confint(lmmatlog, level=89)[2,1]
    lmmatconfintlog89 <- confint(lmmatlog, level=89)[2,2]
    fsestadd <- data.frame(siteslist=sitez[i], cc=unique(fs$cc)[ccstate], meanmat=meanmat, 
                           varmat=varmat, sdmat=sdmat, meanlo=meanlo, varlo=varlo, sdlo=sdlo, meanutah=meanutah, 
                           meangdd=meangdd, matslope=coef(lmmat)["mat"], matslopese=lmmatse, matslopeconfint11=lmmatconfint11, 
                           matslopeconfint89=lmmatconfint89,
                           matslopelog=coef(lmmatlog)["log(mat)"], matslopelogse=lmmatlogse, matslopelogconfint11=lmmatconfintlog11,
                           matslopelogconfint89=lmmatconfintlog89,
                           meanmatlo=meanmatlo,
                           varmatlo=varmatlo, sdmatlo=sdmatlo)
    fsest <- rbind(fsest, fsestadd)
  }
}    

meanhere <- aggregate(fsest[c("meanmat", "varmat", "sdmat", "meanmatlo", "varmatlo", "sdmatlo", "meanlo", "varlo", "sdlo", "meanutah", "meangdd",
                              "matslope", "matslopese", "matslopelog", "matslopelogse")], fsest["cc"], FUN=mean)
sdhere <- aggregate(fsest[c("meanmat", "varmat", "meanmatlo", "varmatlo", "meanlo", "varlo", "meanutah", "meangdd", "matslope")],
                    fsest["cc"], FUN=sd)

#       cc     meanmat   varmat   sdmat  meanmatlo varmatlo   sdmatlo   meanlo    varlo     sdlo meanutah  meangdd
# 1950-1970 7.629376 1.237705 1.1112462  7.810838 1.336977 1.1276187 119.4534 63.49537 7.862192 2093.090 91.08116
# 1990-2010 8.832890 0.901294 0.9488819  7.455157 1.037334 0.9979342 114.0345 31.88875 5.493083 2331.186 79.82799
#   matslope matslopese matslopelog matslopelogse
# -4.651015   1.211960  -0.2883085    0.07817823
# -2.696109   1.182926  -0.2033251    0.09324998


fsest$matslopelog_exp <- exp(fsest$matslopelog)

write.csv(fsest, file="output/fsylestimates_withlog_1950-1970_1990-2010.csv", row.names = FALSE)

## Also get the difference for each site across two time periods
# This is to compare to sims better

fsest.sitediffs <- data.frame(siteslist=numeric(), matdiff=numeric(), matlodiff=numeric(), diffslope=numeric(),
                              varlodiff=numeric(), varlodiffper=numeric(), varmatdiffper=numeric())

for(i in c(1:length(sitez))){ # i <- 1
  subby <- subset(fsest, siteslist==sitez[i])
  precc <- subset(subby, cc=="1950-1960")
  postcc <- subset(subby, cc=="2000-2010")
  matdiff <- precc$meanmat-postcc$meanmat
  matlodiff <- precc$meanmatlo-postcc$meanmatlo
  diffslope <- precc$matslope-postcc$matslope
  diffslopelog <- precc$matslopelog-postcc$matslopelog
  varlodiff <- precc$varlo-postcc$varlo
  varlodiffper <- postcc$varlo/precc$varlo
  varmatdiffper <- postcc$varmat/precc$varmat
  fsest.sitediffs.add <- data.frame(siteslist=sitez[i], matdiff=matdiff,matlodiff=matlodiff, diffslope=diffslope,
                                    diffslopelog,
                                    varlodiff=varlodiff, varlodiffper=varlodiffper, varmatdiffper=varmatdiffper)
  fsest.sitediffs <- rbind(fsest.sitediffs, fsest.sitediffs.add)
}

fsest.sitediffs$daysperC <- fsest.sitediffs$diffslope/fsest.sitediffs$matdiff
fsest.sitediffs$daysperClog <- fsest.sitediffs$diffslopelog/fsest.sitediffs$matdiff


########################
## Stuff in the supp ##
########################

# estimate change in variance
# take the means, then calculates the percent change between pre and post cc 
1-meanhere$varlo[2]/meanhere$varlo[1]
1-mean.pepsims$var.lo.postcc[1]/mean.pepsims$var.lo.precc[1]
1-median.pepsims$var.lo.postcc[1]/median.pepsims$var.lo.precc[1]

# in supp 'we estimated a decline in sensitivity'
mean(fsest.sitediffs$daysperC)
sd(fsest.sitediffs$daysperC)/sqrt(45) # days per C mean SE fs

# in supp 'given X warming'
mean(fsest.sitediffs$matdiff)
mean(fsest.sitediffs$matdiff)/sqrt(45)

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
1 - mean(fsest.sitediffs$varlodiffper)
sd(fsest.sitediffs$varlodiffper)/sqrt(45)
1 - mean.pepsims$varlodiffper[1]
sd.pepsims$varlodiffper[1]/(sqrt(45))

# closer look at variance differences due to when we take the mean ...
# when we take the diff at each site, we exacerbate outliers, this explains (I think) the difference between the two ways we calculate variance
pepsims.1d <- subset(pepsims, degwarm==1)
hist(pepsims.1d$varlodiff)
hist(pepsims.1d$var.lo.precc)
hist(pepsims.1d$var.lo.postcc)

# other stuff ... (not using currently)
mean(fsest.sitediffs$diffslope)
sd(fsest.sitediffs$diffslope)/sqrt(45)
mean(fsest.sitediffs$matlodiff)
sd(fsest.sitediffs$matdiff)
sd(fsest.sitediffs$matlodiff)
sd(fsest.sitediffs$daysperC)

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
realdat.diff <- mean(fsest.sitediffs$diffslope) 
points(abs(mean(fsest.sitediffs$matdiff)), realdat.diff, cex=cexhere, pch=17, col="salmon")
realdatse <- sd(fsest.sitediffs$diffslope)/sqrt(45)
lines(x=rep(abs(mean(fsest.sitediffs$matdiff)), 2), y=c(realdat.diff-realdatse, realdat.diff+realdatse),
      col="salmon")
# par(xpd=TRUE) # so I can plot legend outside
legend("topright", pch=c(17, 19), col=c("salmon", "darkblue"), legend=c("European data (PEP725)", "Simulations with constant cues"),
       cex=1, bty="n")
dev.off()
