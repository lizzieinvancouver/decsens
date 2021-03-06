### Started 1 May 2020 ###
### By Lizzie ###

## Moving some of the extra bits added to decsensSims.R over the months ##
## so I can keep decsensSims.R focused on main bits in paper ##

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

##################################
## Sims to look at correlations ##
##################################

# Step 1: Set up years, days per year, temperatures, required GDD (fstar)
daysperyr <- 60
yearz <- 30
sitez <- 45 # reps
degreez <- c(0, 0.5, 1, 1.5, 2, 2.5, 4, 7)
sigma <- 4
basetemp <- 6
fstar <- 150


# Step 2 -- correlation style!: Build the data and calculate sensitivities with correlations (instead of regression)...
# plus truncated regression for plotting
dfcorr <- data.frame(degwarm=numeric(), rep=numeric(), simplecorr=numeric(), logcorr=numeric(), percorr=numeric(),
    simplecorr.trunc=numeric(), logcorr.trunc=numeric(), simplelm.trunc=numeric(), loglm.trunc=numeric())

for (i in degreez){
   for (j in 1:sitez){
       yearly_expected_temp <- rep(basetemp, yearz)
       daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(daysperyr, basetemp + i, sigma)) 
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > fstar)))
       yearly_temp <- colMeans(daily_temp)
       yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfcorradd <- data.frame(degwarm=i, rep=j, simplecorr=cor(yearly_temp,leafout_date),
           logcorr=cor(log(yearly_temp),log(leafout_date)),
           percorr=cor(per_yearly_temp, per_leafout_date),
           simplecorr.trunc=cor(yearly_temp_trunc,leafout_date),
           logcorr.trunc=cor(log(yearly_temp_trunc),log(leafout_date)),
           simplelm.trunc=coef(lm(leafout_date~yearly_temp_trunc))[2],
           loglm.trunc=coef(lm(log(leafout_date)~log(yearly_temp_trunc)))[2])
       dfcorr <- rbind(dfcorr, dfcorradd)
    }
}

plot(simplecorr~degwarm, data=dfcorr, pch=16, ylab="Sensitivity (days/C or log(days)/log(C)", xlab="Degree warming",
    col="gray", ylim=c(-1, -0.1))
points(simplecorr.trunc~degwarm, data=dfcorr)
points(logcorr~degwarm, data=dfcorr, col="dodgerblue")
points(percorr~degwarm, data=dfcorr, col="firebrick")

plot(abs(simplecorr)~degwarm, data=dfcorr, col="lightgrey",
    ylab="Abs(Sensitivity (days/C or log(days)/log(C))", xlab="Degree warming")
dfcorr$degwarmJitter <- dfcorr$degwarm + 0.05
points(abs(logcorr)~degwarmJitter, data=dfcorr, col="dodgerblue", cex=0.8)


#############################
## Plotting (correlations) ##
############################

# Summarize the sims
mean.sims <- aggregate(dfcorr[c("simplelm.trunc", "loglm.trunc")], dfcorr["degwarm"], FUN=mean)
sd.sims <- aggregate(dfcorr[c("simplelm.trunc", "loglm.trunc")], dfcorr["degwarm"], FUN=sd)

# Get data ...
dfpep <- read.csv("pep_analyses/output/bpenestimates_withlog_1950_2000.csv", header=TRUE)

# Get means and SD
mean.betpen <- aggregate(dfpep[c("meanmat", "lomatslope", "lomatslopelog")], dfpep["cc"], FUN=mean)
sd.betpen <- aggregate(dfpep[c("meanmat", "lomatslope", "lomatslopelog")], dfpep["cc"],  FUN=sd)

tempdiff <- mean.betpen$meanmat[which(mean.betpen$cc=="2000-2010")]-
    mean.betpen$meanmat[which(mean.betpen$cc=="1950-1960")]
tempdiffplot <- c(0, tempdiff)

# PEP plus sims plot with temperature truncated to leafout date version
library(grDevices)
colz <- c("blue4", "violetred4", "blue1", "violetred1")
colzalpha <- adjustcolor(colz, alpha.f = 0.7)

cexhere <- 0.75
cexhere <- 1.2
cextext <- 0.75
jitterpep <- -0.04
pdf(file.path("figures/basicsimsandpep.trunc.pdf"), width = 7.5, height = 5.5)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.25, 2.25), ylim=c(-6.6, 4),yaxt="n",
     ylab=expression(paste("Estimated sensitivity"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="",cex.lab=1.2)
axis(2,seq(-6,0,1),las=2)
# abline(h=0, lty=2, col="darkgrey")
tempsteps <- 5
tempdiffplot <- c(0,1)
for(i in 1:tempsteps){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$simplelm.trunc[i]
  sdhere <- sd.sims$simplelm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[1])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[1])
}
for(i in 1:tempsteps){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$loglm.trunc[i]
  sdhere <- sd.sims$loglm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[2])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[2])
}
for(i in 1:length(unique(mean.betpen$cc))){#i=2
  pos.x <- tempdiffplot[i]+jitterpep
  pos.y <- mean.betpen$lomatslope[i]
  sdhere <- sd.betpen$lomatslope[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[3])
  points(pos.x, pos.y, cex=cexhere, pch=17, col=colzalpha[3])
  text(pos.x + 0.15, pos.y-1.1, labels=unique(mean.betpen$cc)[i], 
       cex=cextext, col=colzalpha[3])
}
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$lomatslopelog[i]
  sdhere <- sd.betpen$lomatslopelog[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[4])
  points(pos.x, pos.y, cex=cexhere, pch=17, col=colzalpha[4])
  text(pos.x + 0.17, pos.y, labels=unique(mean.betpen$cc)[i], cex=cextext, col=colzalpha[4])
}

text(1.1,-6,"Simulations",col="darkgrey")
text(1.1,-6.4,"Observations",col="darkgrey")
text(1.5,-5.5,expression(paste("days / ", degree, "C")),col="darkgrey")
text(2.02,-5.5,expression(paste("log(days) / log(", degree, "C)")),col="darkgrey")

points(1.5, -6, cex=1.7, pch=19, col=colzalpha[1])
points(1.5, -6.4, cex=1.7, pch=17, col=colzalpha[3])
points(2, -6, cex=1.7, pch=19, col=colzalpha[2])
points(2, -6.4, cex=1.7, pch=17, col=colzalpha[4])

lines(c(1.32,1.32),c(-6.6,-5.3),col="darkgrey")
lines(c(1.7,1.7),c(-6.6,-5.3),col="darkgrey")
lines(c(0.85,2.3),c(-5.75,-5.75),col="darkgrey")
lines(c(0.85,2.3),c(-6.2,-6.2),col="darkgrey")

dev.off()


######################################
## Sims to look at imperfect models ##
######################################
# Below sets up 60 days that could lead to leafout, then cuts out the first 20 days when calculating temp
# You can also adjust yearz (as I did to create basicsims.test.10yr.pdf)

# Step 1: Set up years, days per year, temperatures, required GDD (fstar)
daysperyr <- 60
yearz <- 30
sitez <- 45 # reps
degreez.maintext <-  c(0, 0.5, 1, 1.5, 2, 2.5, 4, 7)
sigma <- 4
basetemp <- 6
fstar <- 150

# Step 2: Build the data and calculate sensitivities (note that alpha_1 -- spring temperature increase is set to 0)
df <- data.frame(degwarm=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric(),
    simplelm.trunc=numeric(), loglm.trunc=numeric())

for (i in degreez.maintext){
   for (j in 1:sitez){
       yearly_expected_temp <- rep(basetemp, yearz)
       daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(daysperyr, basetemp + i, sigma)) 
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > fstar)))
       daily_temp <- daily_temp[20:60,]
       yearly_temp <- colMeans(daily_temp)
       yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfadd <- data.frame(degwarm=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
           perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
           simplelm.trunc=coef(lm(leafout_date~yearly_temp_trunc))[2],
           loglm.trunc=coef(lm(log(leafout_date)~log(yearly_temp_trunc)))[2])
       df <- rbind(df, dfadd)
    }
}


# Step 3: Summarize the sims
mean.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["degwarm"], FUN=mean)
sd.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["degwarm"], FUN=sd)

# Step 4: Plot
cexhere <- 0.95
pdf(file.path("figures/basicsims.test.twothirdwindow.pdf"), width = 6, height = 8)
par(mfrow=c(2,1))
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-5, 0.5),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C to leafout)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="",
     bty="l", mgp=c(1.5,.5,0), tck=-.01)
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$simplelm.trunc[i]
  sdhere <- sd.sims$simplelm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$loglm.trunc[i]
  sdhere <- sd.sims$loglm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("salmon", "darkblue"), legend=c("Using logged variables", "Simple linear regression"),
   cex=1, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-5, 0.5),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C over window)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="",
     bty="l", mgp=c(1.5,.5,0), tck=-.01)
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
legend("bottomright", pch=c(19, 19), col=c("salmon", "darkblue"), legend=c("Using logged variables", "Simple linear regression"),
   cex=1, bty="n")
dev.off()
