## Started 28 December 2019 ##
## By Lizzie ##

## Simulations of winter to spring temperatures with fstar (GDD) required increasing when chilling is low ##
## ... and! delayed leafout due to photoperiod cue ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# set.seed(113)

# Setting working directory (only needed to write out figures)
if(length(grep("ailene", getwd()))>0) { 
  setwd("~/Documents/GitHub/decsens/analyses")
} else
  setwd("~/Documents/git/projects/treegarden/decsens/analyses")

# libraries
require(ggplot2)

#############################
## Chilling x forcing sims ##
#############################
# This code requires higher forcing (aka GDD, aka fstar) for leafout when chilling is low #
# It does not vary when GDD accumulates though, it always just starts accumulating on a certain day (daystostart) #

# Step 1: Set up years, days per year, temperatures, required GDD (fstar), required chill (cstar) and how much higher fstar is when cstar is not met
dayswinter <- 120
daysspring <- 90
wintertemp <- 1
springtemp <- 2
springtempinc <- 0.1
sigma <- 4
fstar <- 200
cstar <- 110
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart
yearz <- 50 # used to be 20, adjusted to 50 to compare to Jonathan's code 
sitez <- 45
simsnum <- 30
degreez <- seq(0, 7, length.out=simsnum) # warming -- applied evenly across `winter' and `spring'

if(FALSE){
# Saving some intermediate progress .... setting up the seasonal temps and looking at chill and GDD calculated from them 
df <- data.frame(degwarm=numeric(), rep=numeric(), chill=numeric(), gdd=numeric())

for (i in degreez){
   for (j in 1:sitez){
       yearly_temp <- rep(0, yearz)
       daily_temp <- sapply(yearly_temp, function(x) c(rnorm(dayswinter, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma) + c(1:daysspring)*springtempinc)))
       chill <- daily_temp
       chill[(chill)<0] <- 0
       chill[(chill)>5] <- 0
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       dfadd <- data.frame(degwarm=i, rep=j, chill=mean(colSums(chill[1:dayswinter,], na.rm = TRUE)),
           gdd=mean(colSums(gdd[dayswinter:nrow(gdd),], na.rm = TRUE)))
       df <- rbind(df, dfadd)
    }
}

plot(daily_temp[,1]~c(1:nrow(daily_temp)), type="l")
par(mfrow=c(1,2))
plot(chill~degwarm, df)
abline(h=cstar, col="red", lwd=2)
plot(gdd~degwarm, df)

# Saving some more intermediate progress ... this sets up a varying fstar when chilling is below cstar
gddreq <- c()
leafout_date <- c()
for (k in 1:ncol(chill)){
    chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
    gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
    if (chillsum[dayswinter,k]>cstar) {
        gddreq[k] <- fstar
        } else {
        gddreq[k] <- fstar + (cstar-chillsum[dayswinter,k])*fstaradjforchill
        }
    leafout_date[k] <- min(which(gddsum[,k] > gddreq[k]))
    }

plot(leafout_date~colMeans(daily_temp))
plot(leafout_date~gddreq)
}

## Step 2: Now I put together the seasonal temps, varying fstar (increases when chill is low) and calculate the sensitivities

df <- data.frame(degwarm=numeric(), rep=numeric(), chill=numeric(), fstar=numeric(), simplelm=numeric(),
    loglm=numeric(), perlm=numeric(), propryrschillmet=numeric(), meangddsum=numeric()) 

# yearlytemp <- "postwinter"
yearlytemp <- "alltemps"

plotonesite <- FALSE # this doesn't work when degreez is too big (it worked for length(degreez)=10, but not 30)
if(plotonesite){
# plot the simulated data and simple linear models from one site
site2plot = 11 # arbirarily pick a site to plot! 
figname<-paste("figures/simsiteplots/decsensplot_warm","_site",site2plot,".pdf",sep="")
pdf(figname, width = 12, height = 7)
par(mfrow=c(2, length(degreez)/2))
}
for (i in degreez){
   for (j in 1:sitez){
       yearly_temp <- rep(0, yearz) # set up for sapply
       daily_temp <- sapply(yearly_temp, function(x) c(rnorm(dayswinter, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma) + c(1:daysspring)*springtempinc))) # plot(daily_temp[,1] ~ c(1:length(daily_temp[,1])))
       chill <- daily_temp
       chill[(chill)<0] <- 0
       chill[(chill)>5] <- 0
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       gddreq <- c()
       leafout_date <- c()
       for (k in 1:ncol(chill)){
           chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
           if (chillsum[dayswinter,k]>cstar) {
           gddreq[k] <- fstar
           } else {
           gddreq[k] <- fstar + (cstar-chillsum[dayswinter,k])*fstaradjforchill
           }
           leafout_date[k] <- min(which(gddsum[,k] > gddreq[k])) # leafout_date unit of days *after* dayswinter
           # calculating a few means to include in dataframe
           meanchill <- mean(chillsum[dayswinter,]) 
           meanfstar <- mean(gddreq)
           chillmet<-length(which(chillsum[dayswinter,]>cstar))/yearz
           meangddsum<- mean(gddsum)
           }
           yearly_tempall <- colMeans(daily_temp)
           yearly_temppostwinter <- colMeans(daily_temp[dayswinter:nrow(daily_temp),])
           if(yearlytemp=="postwinter") {
               yearly_temp <- yearly_temppostwinter
               } else {
               yearly_temp <- yearly_tempall
               }
           per_leafout_date <- leafout_date/mean(leafout_date)
           per_yearly_temp <- yearly_temp/mean(yearly_temp)
           # plot(yearly_temp, leafout_date, pch=20)
           dfadd <- data.frame(degwarm=i, rep=j, chill=meanchill, fstar=meanfstar,     
               simplelm=coef(lm(leafout_date~yearly_temp))[2],
               loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
               perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
               propryrschillmet = chillmet,meangddsum= meangddsum)
               if(plotonesite){
                if(j==site2plot){
                  plot(yearly_temp, leafout_date,pch=16,col="darkgreen", bty="l",
                      xlab="Temperature",ylab="Leafout", main = paste(round(i, 2),"warming_site",j, sep=" "))
                  m<-lm(leafout_date~yearly_temp)
                  if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
                }
                }
           df <- rbind(df, dfadd)     
       }
   }
if(plotonesite){
dev.off()
}

if(FALSE){
# Make a more compact version of the above, to share with JA

df <- data.frame(degwarm=numeric(), rep=numeric(), chill=numeric(), fstar=numeric(), simplelm=numeric())
    
for (i in degreez){
   for (j in 1:sitez){
       daily_temp <- sapply(rep(NA, yearz), function(x) c(rnorm(dayswinter, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma)+ c(1:daysspring)*springtempinc)))
       chill <- daily_temp
       chill[(chill)<0] <- 0
       chill[(chill)>5] <- 0
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       gddreq <- c()
       leafout_date <- c()
       for (k in 1:ncol(chill)){
           chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
           if (chillsum[dayswinter,k]>cstar) {
           gddreq[k] <- fstar
           } else {
           gddreq[k] <- fstar + (cstar-chillsum[dayswinter,k])*fstaradjforchill
           }
           leafout_date[k] <- min(which(gddsum[,k] > gddreq[k])) # leafout_date unit of days *after* dayswinter
           meanchill <- mean(chillsum[dayswinter,])
           meanfstar <- mean(gddreq)
           }
           yearly_temp <- colMeans(daily_temp)
           dfadd <- data.frame(degwarm=i, rep=j, chill=meanchill, fstar=meanfstar,     
               simplelm=coef(lm(leafout_date~yearly_temp))[2],
               loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2])
           df <- rbind(df, dfadd)     
       }
   }

}

plot(simplelm~degwarm, data=df, pch=16, ylab="Sensitivity (days/C or log(days)/log(C)", xlab="Degree warming")
points(loglm~degwarm, data=df, col="dodgerblue")
points(perlm~degwarm, data=df, col="firebrick")

plot(propryrschillmet~degwarm, data=df, pch=16, ylab="Number of years (out of 30) when chilling is met", xlab="Degree warming")
plot(fstar~degwarm, data=df, pch=16, ylab="GDD", xlab="Degree warming")
plot(meangddsum~degwarm, data=df, pch=16, ylab="GDD", xlab="Degree warming")



##############
## Plotting ##
##############

mean.sims <- aggregate(df[c("simplelm", "loglm", "perlm","propryrschillmet", "fstar","meangddsum")], df["degwarm"], FUN=mean)
sd.sims <- aggregate(df[c("simplelm", "loglm", "perlm","propryrschillmet", "fstar","meangddsum")], df["degwarm"], FUN=sd)

cexhere <- 0.95
cextext <- 0.8
pdf(file.path("figures/shiftingcuessims_2panels.pdf"), width = 6, height = 8)
par(mfrow=c(2,1),mar=c(5,5,2,5))
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-15, 5),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="", bty="l", mgp=c(1.5,.5,0), tck=-.01)
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$simplelm[i]
  sdhere <- sd.sims$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$loglm[i]
  sdhere <- sd.sims$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
# par(xpd=TRUE) # so I can plot legend outside
legend("bottomright", pch=c(19, 19), col=c( "salmon","darkblue"), legend=c("Using logged variables","Simple linear regression"),
   cex=cextext, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-0.1, 1),mgp=c(1.5,.5,0), tck=-.01,xaxs="i",yaxs = "i",
     ylab="Proportion years when chilling is met",
     xlab=expression(paste("Warming (", degree, "C)")), bty="u",main="")
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$propryrschillmet[i]
  sdhere <- sd.sims$propryrschillmet[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkgray")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkgray")
}
par(new = TRUE)

plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(200,300),yaxt="n", ylab="",xaxt="n", xlab="", bty="u",mgp=c(1.5,.5,0), tck=-.01)
axis(side = 4,mgp=c(1.5,.5,0), tck=-.01)
mtext(expression(paste("Thermal sum required for leafout (", degree, "C)"), sep=""),side=4, adj=.5, line=2)

for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$fstar[i]
  sdhere <- sd.sims$fstar[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkred")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkred")
}
text(mean.sims$degwarm[1] + 1.6, mean.sims$fstar[1] + 6,
     labels=expression(paste("Thermal sum (", degree, "C)"), sep=""), cex=cextext,
      col="darkred")
text(mean.sims$degwarm[1] + 0.8, mean.sims$fstar[1] + 80,
     labels=expression(paste("Chilling met"), sep=""), cex=cextext,
      col="darkgray")

dev.off()



################################################
## From Jonathan, understanding shifting cues ##
################################################
# taken from auerbach_sep2020.R

dayswinter <- 120
daysspring <- 90
wintertemp <- 1
springtemp <- 2
springtempinc <- 0.1
sigma <- 4
fstar <- 200
cstar <- 110
yearz <- 50
sitez <- 45
simsnum <- 11
degreez <- 0:10

df <- data.frame(degwarm=numeric(), rep=numeric(), chill=numeric(),
    fstar=numeric(), simplelm=numeric(), fstaradjforchill  = numeric(), gddreq = numeric())

for(fstaradjforchill in c(0,1,3,9)) {
for (i in degreez){
    for (j in 1:sitez){
        daily_temp <- sapply(rep(NA, yearz), function(x)
          c(rnorm(dayswinter, wintertemp + i, sigma),
            (rnorm(daysspring, springtemp + i , sigma)+
               c(1:daysspring)*springtempinc)))
       
        chill <- daily_temp
        chill[(chill)<0] <- 0
        chill[(chill)>5] <- 0
        gdd <- daily_temp
        gdd[(gdd)<0] <- 0
        gddreq <- c()
        leafout_date <- c()
        for (k in 1:ncol(chill)){
            chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
            gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
         
          if (chillsum[dayswinter,k]>cstar) {
            gddreq[k] <- fstar
            } else {
            gddreq[k] <- fstar + (cstar-chillsum[dayswinter,k])*fstaradjforchill
            }
           
            leafout_date[k] <- min(which(gddsum[,k] > gddreq[k]))

            meanchill <- mean(chillsum[dayswinter,])
            meanfstar <- mean(gddreq)
            }
            yearly_temp <- colMeans(daily_temp)
            if(FALSE){ # this can cause some negative yearly_temps, and then you can't take the log
            for(col in 1:length(yearly_temp)) yearly_temp[col] <-
                mean(daily_temp[1:leafout_date[col],col])
            }
            dfadd <- data.frame(degwarm=i,
                                rep=j,
                                chill=meanchill,
                                fstar=meanfstar,
                                simplelm=coef(lm(leafout_date~yearly_temp))[2],
                                loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
                                fstaradjforchill = fstaradjforchill,
                                gddreq = mean(gddreq))
            df <- rbind(df, dfadd)
        }
    }
}

ggplot(df) +
  theme_bw() +
  aes(factor(degwarm), loglm) +
  geom_boxplot() +
  labs(x = "temperature increase", y = "sensitivity (log)",
       title = "Fig. 1") +
  facet_wrap(~ fstaradjforchill, scales = "free",
             labeller = labeller(.cols = label_both))

ggplot(df) +
  theme_bw() +
  aes(factor(degwarm), simplelm) +
  geom_boxplot() +
  labs(x = "temperature increase", y = "sensitivity",
       title = "Fig. 2") +
  facet_wrap(~ fstaradjforchill, scales = "free",
             labeller = labeller(.cols = label_both))

ggplot(df) +
  theme_bw() +
  aes(factor(degwarm), gddreq) +
  geom_boxplot() +
  labs(x = "temperature increase", y = "average required gdd",
       title = "Fig. 3") +
  facet_wrap(~ fstaradjforchill, scales = "free",
             labeller = labeller(.cols = label_both))
 
ggplot(df) +
  theme_bw() +
  aes(gddreq, loglm) +
  geom_point() +
  geom_smooth() +
  labs(x = "average required gdd", y = "sensitivity (log)",
       title = "Fig. 4") +
  facet_wrap(~ fstaradjforchill, scales = "free",
             labeller = labeller(.cols = label_both))

ggplot(df) +
  theme_bw() +
  aes(gddreq, simplelm) +
  geom_point() +
  geom_smooth() +
  labs(x = "average required gdd", y = "sensitivity",
       title = "Fig. 5") +
  facet_wrap(~ fstaradjforchill, scales = "free",
             labeller = labeller(.cols = label_both))



##########################
## Non-identifiability  ##
##########################

# Simplify code to look for years with same DOY for different cstar and fstar values

dayswinter <- 120
daysspring <- 90
wintertemp <- 1
springtemp <- 2
springtempinc <- 0.1
sigma <- 4 
fstar <- 200
cstar <- 110
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart
yearz <- 2000
simsnum <- 1
degreez <- seq(0, 0, length.out=simsnum) 

df.noninf <- data.frame(degwarm=numeric(), year=numeric(), chill=numeric(), fstar=numeric(), 
    meantemp=numeric(), meantempspr=numeric(), leafoutdoy=numeric()) 

for (i in degreez){
       yearly_temp <- rep(0, yearz)
       daily_temp <- sapply(yearly_temp, function(x) c(rnorm(dayswinter, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma)+ c(1:daysspring)*springtempinc)))
       chill <- daily_temp
       chill[(chill)<0] <- 0
       chill[(chill)>5] <- 0
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       gddreq <- c()
       leafout_date <- c()
       for (k in 1:ncol(chill)){
           chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
           if (chillsum[dayswinter,k]>cstar) {
           gddreq[k] <- fstar
           } else {
           gddreq[k] <- fstar + (cstar-chillsum[dayswinter,k])*fstaradjforchill
           }
           leafout_date[k] <- min(which(gddsum[,k] > gddreq[k])) # leafout_date unit of days *after* dayswinter
           dfadd <- data.frame(degwarm=i, year=k, chill=chillsum[dayswinter,k], fstar=gddreq[k],
               meantemp=mean(daily_temp[,k]), meantempspr=mean(daily_temp[dayswinter:nrow(daily_temp),k]),
               leafoutdoy=leafout_date[k])
           df.noninf <- rbind(df.noninf, dfadd)     

           }
       }


ggplot(data=df.noninf, aes(x=fstar, y=leafoutdoy, group=degwarm, color=as.factor(degwarm))) +
   geom_point()

day50 <- subset(df.noninf, leafoutdoy==50)
plot(chill~fstar, data=day50)

# plot by fstar
fstarhere <- fstar
# saved as shiftingcuessims_chillexplore.pdf
par(mfrow=c(1,2))
plot(leafoutdoy ~ meantemp, data=df.noninf, col="gray", main="winter-spring temps")
abline(lm(leafoutdoy ~ meantemp, data=subset(df.noninf, fstar==fstarhere)), col="navy", lwd=2)
abline(lm(leafoutdoy ~ meantemp, data=subset(df.noninf, fstar>fstarhere)), col="darkred", lwd=2)
abline(lm(leafoutdoy ~ meantemp, data=df.noninf), lwd=2)
legend("topright", lty=rep(1,3), col=c("black", "navy", "darkred"), legend=c("all sims", "chill met", "chill not met"),
   cex=1, bty="n")

summary(lm(leafoutdoy ~ meantemp, data=df.noninf))
summary(lm(leafoutdoy ~ meantemp, data=subset(df.noninf, fstar>fstarhere)))
summary(lm(leafoutdoy ~ meantemp, data=subset(df.noninf, fstar==fstarhere)))

plot(leafoutdoy ~ meantempspr, data=df.noninf, col="gray", main="spring temps")
abline(lm(leafoutdoy ~ meantempspr, data=subset(df.noninf, fstar>fstarhere)), col="darkred", lwd=2)
abline(lm(leafoutdoy ~ meantempspr, data=subset(df.noninf, fstar==fstarhere)), col="navy", lwd=2)
abline(lm(leafoutdoy ~ meantempspr, data=df.noninf), lwd=2)

summary(lm(leafoutdoy ~ meantempspr, data=df.noninf))
summary(lm(leafoutdoy ~ meantempspr, data=subset(df.noninf, fstar>fstarhere)))
summary(lm(leafoutdoy ~ meantempspr, data=subset(df.noninf, fstar==fstarhere)))
# So, if chilling is not met and you're using simple lm you should expect a HIGHER sensitivity, not lower...
# and this holds across two temperature windows.


##############################
## Daylength x forcing sims ##
##############################

# There are two scenarios I can imagine where photoperiod cues kick in with climate change
# (1) Plants reach threshold GDD too soon and thus wait for a certain photoperiod
# (2) Chilling is not met, so required GDD goes way up and plants never reach that, so they just leaf out at a certain photoperiod
# I model a simplified version of this below (effectively (1), but it's pretty similar to (2))
# Here I modeled a simple threshold photoperiod (I did this to make it more different than the chill scenarios), but ...
# could easily make it interactive similar to chilling code above.

# Step 1: Set up years, days per year, temperatures, required GDD (fstar), required photoperiod (pstar)

library(geosphere) # for daylengths

dayswinter <- 60 # sims are set up as though starting on Jan 1 (so 60 here means 'winter temps' end in early March)
daysspring <- 90
wintertemp <- 0
springtemp <- 4
springtempinc <- 0.1
photoper <- daylength(45, 1:(dayswinter+daysspring)) # latitude=45
sigma <- 4
fstar <- 200
pstar <- 12 # threshold to leafout
pstarday <- min(which(photoper > pstar))
yearz <- 50
sitez <- 45
simsnum <- 30
degreez <- round(seq(0, 7, length.out=simsnum), digits=1) # warming -- applied evenly across `winter' and `spring'

# Step 2: Run the simulations
dfphoto <- data.frame(degwarm=numeric(), rep=numeric(), meantemp=numeric(), leafoutdoy=numeric(), gddmetday=numeric(), 
    propyrsphoto=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric()) 

for (i in degreez){
   for (j in 1:sitez){
       yearly_temp <- rep(0, yearz) 
       daily_temp <- sapply(yearly_temp, function(x) c(rnorm(dayswinter, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma)+ c(1:daysspring)*springtempinc)))
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       leafout_date <- c()
       gddmetday <- c()
       for (k in 1:ncol(daily_temp)){
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
           gddmetday[k] <- min(which(gddsum[,k] > fstar))
           if (photoper[(gddmetday[k] + dayswinter)] > pstar) { # check if the photoperiod threshold is met by the GDD-met day
            leafout_date[k] <- gddmetday[k]
           } else {
            leafout_date[k] <- pstarday-dayswinter # keep on same day scale as gdd
           }
           }
           yearly_temp <- colMeans(daily_temp)
           per_leafout_date <- leafout_date/mean(leafout_date)
           per_yearly_temp <- yearly_temp/mean(yearly_temp)
           driverbyyear <- NA
           driverbyyear[which(leafout_date==gddmetday)] <- "forcing"
           driverbyyear[which(leafout_date!=gddmetday)] <- "photo"
           driverbyyear[which(leafout_date==(pstarday-dayswinter))] <- "forcing/photo"
           photodriver <- length(which(driverbyyear!="forcing"))/yearz # here I pick only definitive forcing years!
           dfadd <- data.frame(degwarm=i, rep=j, meantemp=mean(yearly_temp),
               leafoutdoy=mean(leafout_date), gddmetday= mean(gddmetday),
               propyrsphoto=photodriver,
               simplelm=coef(lm(leafout_date~yearly_temp))[2],
               loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
               perlm=coef(lm(per_leafout_date~per_yearly_temp))[2])
           dfphoto <- rbind(dfphoto, dfadd)     
       }
   }

# quick look at leafout days
par(mfrow=c(1,1))
pstarday-dayswinter
hist(dfphoto$leafoutdoy)


# plot by whether leafout was due to forcing and/or photoperiod (averaged over years)
# propyrsphoto for better resolution, but even then it is not totally accurate when the gdd met day is same as photo threhold!
dfphoto$leafoutdriver <- NA
dfphoto$leafoutdriver[which(dfphoto$propyrsphoto==0)] <- "all forcing"
dfphoto$leafoutdriver[which(dfphoto$propyrsphoto==1)] <- "all photo" # in current sims does not always occur
dfphoto$leafoutdriver[which(dfphoto$propyrsphoto>0 & dfphoto$propyrsphoto<1)] <- "photo/forcing"
dfphoto$degwarmtext <- paste("warming:", as.factor(dfphoto$degwarm), "C")

if(FALSE){
# saved as shiftingcuessims_photo.pdf
ggplot(data=dfphoto, aes(x=meantemp, y=leafoutdoy, group=leafoutdriver, color=leafoutdriver)) +
   geom_point() +
   geom_smooth(method = "lm", linetype = 2, lwd=0.5, se = FALSE) +
   facet_wrap(.~degwarmtext, scales="free") +
        ylab("Mean day of year") +
    xlab(expression(paste("Mean daily temperature (", degree, "C)"), sep=""))

# saved as shiftingcuessims_photogradient.pdf
ggplot(data=dfphoto, aes(x=meantemp, y=leafoutdoy, color=propyrsphoto)) +
   geom_point() +
   geom_smooth(method = "lm", linetype = 2, lwd=0.5, se = FALSE) +
   facet_wrap(.~degwarmtext, scales="free") +
        ylab("Mean day of year") +
    xlab(expression(paste("Mean daily temperature (", degree, "C)"), sep=""))
}

# plot the means
mean.simsphoto <- aggregate(dfphoto[c("simplelm", "loglm", "perlm","leafoutdoy", "gddmetday", "propyrsphoto")], dfphoto["degwarm"], FUN=mean)
sd.simsphoto <- aggregate(dfphoto[c("simplelm", "loglm", "perlm","leafoutdoy", "gddmetday", "propyrsphoto")], dfphoto["degwarm"], FUN=sd)

cexhere <- 0.95
cextext <- 0.8

pdf(file.path("figures/shiftingcuessims_photo2panel.pdf"), width = 6, height = 8)
par(mfrow=c(2,1), mar=c(5,5,2,5))
plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(-8, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="", bty="l", mgp=c(1.5,.5,0), tck=-.01)
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$simplelm[i]
  sdhere <- sd.simsphoto$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$loglm[i]
  sdhere <- sd.simsphoto$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c( "salmon","darkblue"), legend=c("Using logged variables","Simple linear regression"),
   cex=cextext, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez)+ 0.5)), ylim=c(-0.1, (max(mean.simsphoto$propyrsphoto)+0.1)),
     mgp=c(1.5,.5,0), tck=-.01,xaxs="i",yaxs = "i",
     ylab="Proportion years photoperiod drives leafout",
     xlab=expression(paste("Warming (", degree, "C)")), bty="l",main="")
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$propyrsphoto[i]
  sdhere <- sd.simsphoto$propyrsphoto[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkgray")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkgray")
}
dev.off()



if(FALSE){
# For just one site ...
dfphoto_onesite <- data.frame(degwarm=numeric(), rep=numeric(), yearly_temp=numeric(), leafoutdoy=numeric(),
    gddmetday=numeric(), driver=numeric()) 

for (i in degreez){
       yearly_temp <- rep(0, yearz) 
       daily_temp <- sapply(yearly_temp, function(x) c(rnorm(dayswinter, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma)+ c(1:daysspring)*springtempinc)))
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       leafout_date <- c()
       gddmetday <- c()
       for (k in 1:ncol(daily_temp)){
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[dayswinter:nrow(gdd),x])))
           gddmetday[k] <- min(which(gddsum[,k] > fstar))
           if (photoper[(gddmetday[k] + dayswinter)] > pstar) { # check if the photoperiod threshold is met by the GDD-met day
            leafout_date[k] <- gddmetday[k]
           } else {
            leafout_date[k] <- pstarday-dayswinter # keep on same day scale as gdd
           }
           }
           yearly_temp <- colMeans(daily_temp)
           per_leafout_date <- leafout_date/mean(leafout_date)
           per_yearly_temp <- yearly_temp/mean(yearly_temp)
           driverbyyear <- NA
           driverbyyear[which(leafout_date==gddmetday)] <- "forcing"
           driverbyyear[which(leafout_date!=gddmetday)] <- "photo"
           driverbyyear[which(gddmetday==(pstarday-dayswinter))] <- "forcing/photo"
           dfadd <- data.frame(degwarm=rep(i, yearz), yearly_temp=yearly_temp,
               leafoutdoy=leafout_date, gddmetday=gddmetday,
               driver=driverbyyear)
           dfphoto_onesite <- rbind(dfphoto_onesite, dfadd)     
       }


dfphoto_onesite$degwarmtext <- paste("warming:", as.factor(dfphoto_onesite$degwarm), "C")

# saved as shiftingcuessims_photo_onesite.pdf
ggplot(data=dfphoto_onesite, aes(x=yearly_temp, y=leafoutdoy, group=driver, color=driver)) +
   geom_point() +
   geom_smooth(method = "lm", linetype = 2, lwd=0.5, se = FALSE) +
   facet_wrap(.~degwarmtext, scales="free") +
        ylab("Day of year") +
    xlab(expression(paste("Mean daily temperature (", degree, "C)"), sep=""))
}

###################################
## Plots chilling and daylength  ##
##          sims together        ##
###################################


pdf(file.path("figures/shiftingcuessims_4panels.pdf"), width = 12, height = 8)
par(mfrow=c(2,2), mar=c(5,5,2,5))
plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(-10, 5),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="", bty="l", mgp=c(1.5,.5,0), tck=-.01)
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$simplelm[i]
  sdhere <- sd.sims$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$loglm[i]
  sdhere <- sd.sims$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("salmon","darkblue"), legend=c("Using logged variables","Simple linear regression"),
   cex=cexhere, bty="n")


plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(-6, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="", bty="l", mgp=c(1.5,.5,0), tck=-.01)
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$simplelm[i]
  sdhere <- sd.simsphoto$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$loglm[i]
  sdhere <- sd.simsphoto$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c( "salmon","darkblue"), legend=c("Using logged variables","Simple linear regression"),
   cex=cexhere, bty="n")

plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(-0.1, 1),mgp=c(1.5,.5,0), tck=-.01,xaxs="i",yaxs = "i",
     ylab="Proportion years when chilling is met",
     xlab=expression(paste("Warming (", degree, "C)")), bty="u",main="")
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$propryrschillmet[i]
  sdhere <- sd.sims$propryrschillmet[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkgray")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkgray")
}
par(new = TRUE)
plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez) + 0.5)), ylim=c(200,300),yaxt="n", ylab="",xaxt="n", xlab="", bty="u",mgp=c(1.5,.5,0), tck=-.01)
axis(side = 4,mgp=c(1.5,.5,0), tck=-.01)
mtext(expression(paste("Thermal sum required for leafout (", degree, "C)"), sep=""), side=4, adj=.5, line=2, cex=cexhere)

for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$fstar[i]
  sdhere <- sd.sims$fstar[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkred")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkred")
}
text(mean.sims$degwarm[1] + 1.6, mean.sims$fstar[1] + 6,
     labels=expression(paste("Thermal sum (", degree, "C)"), sep=""), cex=cextext,
      col="darkred")
text(mean.sims$degwarm[1] + 0.8, mean.sims$fstar[1] + 80,
     labels=expression(paste("Chilling met"), sep=""), cex=cextext,
     col="darkgray")


plot(x=NULL,y=NULL, xlim=c(-0.5, (max(degreez)+ 0.5)), ylim=c(-0.1, (max(mean.simsphoto$propyrsphoto)+0.1)),
     mgp=c(1.5,.5,0), tck=-.01,xaxs="i",yaxs = "i",
     ylab="Proportion years photoperiod drives leafout",
     xlab=expression(paste("Warming (", degree, "C)")), bty="l",main="")
for(i in 1:length(unique(mean.simsphoto$degwarm))){
  pos.x <- mean.simsphoto$degwarm[i]
  pos.y <- mean.simsphoto$propyrsphoto[i]
  sdhere <- sd.simsphoto$propyrsphoto[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkgray")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkgray")
}
dev.off()
