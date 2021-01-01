## Started 31 December 2020 ##
## By Lizzie ##

## Making year-long climate data to play around with sliding window... ## 
## Created in part from decsensSimsMultCues.R ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/treegarden/decsens/analyses/pep_analyses") 
} else setwd("~/Documents/git/decsens/analyses/pep_analyses")

plottinghere=FALSE # make some plots of the data?

dayswinter <- 100
daysspring <- 90
dayssummer <- 85
daysfall <- 90
dayz <- sum(dayswinter, daysspring, dayssummer, daysfall)
wintertemp <- 1
springtemp <- 2
springtempinc <- 0.1
summertemp <- springtemp+springtempinc*daysspring+1
sigma <- 4
fstar <- 200
cstar <- 110
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart
yearz <- 30 
sitez <- 45
simsnum <- 1
degreez <- seq(0, 7, length.out=simsnum) # warming -- applied evenly across `winter' and `spring'

dfbb <- data.frame(site=numeric(), year=numeric(), gddreq=numeric(),
    leafout_date=numeric())
dfclim <- data.frame(site=numeric(), year=numeric(), day=numeric(), daily_temp=numeric())

# Climate data goes dayswinter/2, daysspring, summer, dayssummer, daysfall, dayswinter/2
# Leafout_date is reported in df as days since start of spring PLUS dayswinter/2
for (i in degreez){
   for (j in 1:sitez){
       daily_temp <- sapply(rep(NA, yearz), function(x) c(rnorm(dayswinter/2, wintertemp + i, sigma),
           (rnorm(daysspring, springtemp + i , sigma) + c(1:daysspring)*springtempinc),
           rnorm(dayssummer, summertemp + i, sigma),
           (rnorm(daysfall, summertemp + i , sigma) - c(1:daysfall)*springtempinc),
           rnorm(dayswinter/2, wintertemp + i, sigma)))
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
           dfclimadd <- data.frame(site=rep(j, dayz), year=rep(k, dayz), day=1:dayz, daily_temp=daily_temp[,k])
           dfclim <- rbind(dfclim, dfclimadd)
           }
           yearly_temp <- colMeans(daily_temp)
           dfbbadd <- data.frame(site=rep(j, yearz), year=1:yearz, gddreq=gddreq, leafout_date=leafout_date+dayswinter/2)
           dfbb <- rbind(dfbb, dfbbadd)

       }
   }

if(plottinghere){
plot(daily_temp[,1]~c(1:dayz), type="l")
hist(dfbb$leafout_date)

require(ggplot2)

ggplot(dfclim, aes(x=as.numeric(day), y=daily_temp)) +
    geom_line() + 
    facet_wrap(.~as.factor(site)) +
    xlab("Day of year") +
    ylab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()
}

if((i==0) && (fstaradjforchill>0)){
write.csv(dfclim, "simmonds_slidingwin/output/fakeclim_wchill_0deg.csv")
write.csv(dfbb, "simmonds_slidingwin/output/fakebb_wchill_0deg.csv")
}
if((i==2) && (fstaradjforchill>0)){
write.csv(dfclim, "simmonds_slidingwin/output/fakeclim_wchill_2deg.csv")
write.csv(dfbb, "simmonds_slidingwin/output/fakebb_wchill_2deg.csv")
}

if((i==0) && (fstaradjforchill==0)){
write.csv(dfclim, "simmonds_slidingwin/output/fakeclim_nochill_0deg.csv")
write.csv(dfbb, "simmonds_slidingwin/output/fakebb_nochill_0deg.csv")
}
if((i==2) && (fstaradjforchill==0)){
write.csv(dfclim, "simmonds_slidingwin/output/fakeclim_nochill_2deg.csv")
write.csv(dfbb, "simmonds_slidingwin/output/fakebb_nochill_2deg.csv")
}

