## Started 7 September 2020 ##
## By Lizzie, based off decsensSims.R ##

## Simulation of the declining sensitivities problem ##
## What happens with no warming where only fstar varies? ##

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
  setwd("~/GitHub/decsens") 
} else if(length(grep("catchamberlain", getwd()))>0) { 
  setwd("~/Documents/git/ospree/analyses/bb_analysis/pep_sims")
} else setwd("~/Documents/git/projects/treegarden/decsens/analyses")

##########################
# The below sets up data #
##########################

# Make some data ... note that this runs simsnum times for 45 sites, via a loop
# which is to say, it's slow!

# Step 1: Set up years, days per year, temperatures, required GDD (fstar)
daysperyr <- 60
yearz <- 30
sitez <- 45 # reps
simsnum <- 40
fstarsims <- seq(100, 300, length.out=simsnum)
sigma <- 4
basetemp <- 4 # alpha_0
dailytempchange <- 0.1 # alpha_1

# Step 2: Build the data and calculate sensitivities 
df <- data.frame(fstar=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric(),
    simplelm.trunc=numeric(), loglm.trunc=numeric())

for (i in fstarsims){
   for (j in 1:sitez){
       daily_temp <- sapply(rep(NA, yearz), function(x) rnorm(daysperyr, basetemp, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange 
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > i)))
       yearly_temp <- colMeans(daily_temp)
       yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfadd <- data.frame(fstar=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
           perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
           simplelm.trunc=coef(lm(leafout_date~yearly_temp_trunc))[2],
           loglm.trunc=coef(lm(log(leafout_date)~log(yearly_temp_trunc)))[2])
       df <- rbind(df, dfadd)
    }
}

dfsm <- data.frame(fstar=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric())

for (i in fstarsims){
   for (j in 1:sitez){
       daily_temp <- sapply(rep(NA, yearz), function(x) rnorm(daysperyr, basetemp, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange 
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > i)))
       yearly_temp <- colMeans(daily_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfadd <- data.frame(fstar=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2])
       dfsm <- rbind(dfsm, dfadd)
    }
}

mean.sims.sm <- aggregate(dfsm[c("simplelm", "loglm")], dfsm["fstar"], FUN=mean)


##############
## Plotting ##
##############

# Summarize the fstar sims
mean.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["fstar"], FUN=mean)
sd.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["fstar"], FUN=sd)

colz <- c("blue4", "violetred4")
colzalpha <- adjustcolor(colz, alpha.f = 0.7)
cexhere <- 0.75
cexhere <- 1.2
cextext <- 0.75
jitterpep <- -0.04
pdf(file.path("figures/basicsims_fstaronly.pdf"), width = 7.5, height = 5.5)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(90, 310), ylim=c(-6.6, 1), yaxt="n",
     ylab=expression(paste("Estimated sensitivity"), sep=""),
     xlab=expression(paste("Thermal sum required (", degree, "C)")), main="", cex.lab=1.2,
     bty="l", mgp=c(1.5,.5,0), tck=-.01)
axis(2,seq(-6,0,1),las=2)
tempsteps <- simsnum
tempdiffplot <- c(0,1)
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$simplelm[i]
  sdhere <- sd.sims$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[1])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[1])
}
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$loglm[i]
  sdhere <- sd.sims$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[2])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[2])
}
legend("topright", pch=c(19, 19), col=c(colzalpha[1], colzalpha[2]), legend=c("linear", "non-linear (logged)"),
   cex=1, bty="n")
dev.off()
