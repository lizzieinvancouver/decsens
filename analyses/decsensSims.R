## Started 4 January 2020 ##
## By Lizzie, see also decsensSimsAuerbach.R and pepvarsim.R ##

## Simulation of the declining sensitivities problem ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

set.seed(113)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("ailene", getwd()))>0) { 
  setwd("~/Documents/GitHub/ospree/analyses/decsens")
} else if(length(grep("lizzie", getwd()))>0) {
  setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/decsens")
} else if(length(grep("Ignacio", getwd()))>0) { 
  setwd("~/GitHub/ospree/analyses/decsens") 
} else if(length(grep("catchamberlain", getwd()))>0) { 
  setwd("~/Documents/git/ospree/analyses/bb_analysis/pep_sims")
} else setwd("~/Documents/git/projects/treegarden/decsens/analyses")

# Make some data ...

# Step 1: Set up years, days per year, temperatures, sampling frequency, required GDD (fstar)
daysperyr <- 60
yearz <- 30
sitez <- 45 # reps
degreez <- c(0, 0.5, 1, 1.5, 2, 2.5, 4, 7)
sigma <- 4
basetemp <- 6
fstar <- 150

# Step 2: Build the data and calculate sensitivities
df <- data.frame(degwarm=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric())

for (i in degreez){
   for (j in 1:sitez){
       yearly_expected_temp <- rep(6, yearz)
       daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(daysperyr, basetemp + i, sigma)) 
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > fstar)))
       yearly_temp <- colMeans(daily_temp)
       # yearly_temp <- rnorm(length(yearly_temp), yearly_temp, 1) # add noise to reduce slopes
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       # yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       dfadd <- data.frame(degwarm=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
           perlm=coef(lm(per_leafout_date~per_yearly_temp))[2])
       df <- rbind(df, dfadd)
    }
}

plot(simplelm~degwarm, data=df, pch=16, ylab="Sensitivity (days/C or log(days)/log(C)", xlab="Degree warming")
points(loglm~degwarm, data=df, col="dodgerblue")
plot(perlm~degwarm, data=df, col="firebrick")


plot(abs(simplelm)~degwarm, data=df, col="lightgrey",
    ylab="Abs(Sensitivity (days/C or log(days)/log(C))", xlab="Degree warming")
df$degwarmJitter <- df$degwarm + 0.05
points(abs(loglm)~degwarmJitter, data=df, col="dodgerblue", cex=0.8)


##############
## Plotting ##
##############

mean.sims <- aggregate(df[c("simplelm", "loglm", "perlm")], df["degwarm"], FUN=mean)
sd.sims <- aggregate(df[c("simplelm", "loglm", "perlm")], df["degwarm"], FUN=sd)


cexhere <- 0.95
pdf(file.path("figures/basicsims.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-6, -0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
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
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Simple linear regression", "Using logged variables"),
   cex=1, bty="n")
dev.off()



########################
## Plotting, plus PEP ##
########################

# Get data ...
dfpep <- read.csv("pep_analyses/output/bpenestimates_withlog_1950_2000.csv", header=TRUE) 

# Get means and SD
mean.betpen <- aggregate(dfpep[c("matslope", "matslopelog", "meanmat")], dfpep["cc"], FUN=mean)
sd.betpen <- aggregate(dfpep[c("matslope", "matslopelog", "meanmat")], dfpep["cc"],  FUN=sd)

tempdiff <- mean.betpen$meanmat[which(mean.betpen$cc=="2000-2010")]-
    mean.betpen$meanmat[which(mean.betpen$cc=="1950-1960")]
tempdiffplot <- c(0, tempdiff)

library(grDevices)
colz <- c("blue4", "violetred4", "blue1", "violetred1")
colzalpha <- adjustcolor(colz, alpha.f = 0.5)

cexhere <- 0.75
pdf(file.path("figures/basicsimsandpep.pdf"), width = 5, height = 3.75)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.25, 2.25), ylim=c(-6, -0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
tempsteps <- 5
for(i in 1:tempsteps){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$simplelm[i]
  sdhere <- sd.sims$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[1])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[1])
  }
for(i in 1:tempsteps){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$loglm[i]
  sdhere <- sd.sims$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[2])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[2])
  }
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslope[i]
  sdhere <- sd.betpen$matslope[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[3])
  points(pos.x, pos.y, cex=cexhere, pch=17, col=colzalpha[3])
  }
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslopelog[i]
  sdhere <- sd.betpen$matslopelog[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[4])
  points(pos.x, pos.y, cex=cexhere, pch=17, col=colzalpha[4])
  }
# par(xpd=TRUE) # so I can plot legend outside
legend("bottomright", pch=c(19, 17, 19, 17), col=colzalpha[c(1,3,2,4)],
       legend=c("Sims: untransformed", "PEP725: untransformed",
       "Sims: logged", "PEP725: logged"), cex=0.75, bty="n") # xjust=1 is not working, arghh
dev.off()
