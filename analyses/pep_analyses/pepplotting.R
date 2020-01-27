## Started 18 January 2020 ##
## By Lizzie ##

## PEP data ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Setting working directory. 
setwd("~/Documents/git/projects/treegarden/decsens/analyses")

# Get data ...
df <- read.csv("pep_analyses/output/bpenestimates_withlog.csv", header=TRUE)
dfswa <- read.csv("pep_analyses/output/swaestimates_withlog.csv", header=TRUE)
df90s <- read.csv("pep_analyses/output/bpenestimates_withlog_1950_1990_2000.csv", header=TRUE)

fs90s <- read.csv("pep_analyses/output/fsylestimates_withlog_1950_1990_2000.csv", header=TRUE)

##############
## Plotting ##
##############

mean.betpen <- aggregate(df[c("matslope", "matslopelog", "meanmat")], df["cc"], FUN=mean)
sd.betpen <- aggregate(df[c("matslope", "matslopelog", "meanmat")], df["cc"],  FUN=sd)

mean.betpen90s <- aggregate(df90s[c("matslope", "matslopelog", "meanmat")], df90s["cc"], FUN=mean)
sd.betpen90s <- aggregate(df90s[c("matslope", "matslopelog", "meanmat")], df90s["cc"],  FUN=sd)

tempdiff <- mean.betpen$meanmat[which(mean.betpen$cc=="2000-2010")]-
    mean.betpen$meanmat[which(mean.betpen$cc=="1950-1960")]
tempdiff90s <- mean.betpen90s$meanmat[which(mean.betpen90s$cc=="1990-2000")]-
    mean.betpen90s$meanmat[which(mean.betpen90s$cc=="1950-1960")]
tempdiffplot <- c(0, tempdiff)
tempdiffplot.w90s <- c(0, tempdiff90s, tempdiff)


cexhere <- 0.95
pdf(file.path("figures/basicpep.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.5, 3), ylim=c(-6, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslope[i]
  sdhere <- sd.betpen$matslope[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslopelog[i]
  sdhere <- sd.betpen$matslopelog[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Simple linear regression", "Using logged variables"),
   cex=1, bty="n")
dev.off()

# With the 90s addedd
cexhere <- 0.95
pdf(file.path("figures/pepbetpenw90s.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.5, 3), ylim=c(-6, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen90s$cc))){
  pos.x <- tempdiffplot.w90s[i]
  pos.y <- mean.betpen90s$matslope[i]
  sdhere <- sd.betpen90s$matslope[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.betpen90s$cc))){
  pos.x <- tempdiffplot.w90s[i]
  pos.y <- mean.betpen90s$matslopelog[i]
  sdhere <- sd.betpen90s$matslopelog[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
# par(xpd=TRUE) # so I can plot legend outside
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Simple linear regression", "Using logged variables"),
   cex=1, bty="n")
dev.off()

## For Fagus sylvatica
mean.fagsyl90s <- aggregate(fs90s[c("matslope", "matslopelog", "meanmat")], fs90s["cc"], FUN=mean)
sd.fagsyl90s <- aggregate(fs90s[c("matslope", "matslopelog", "meanmat")], fs90s["cc"],  FUN=sd)

tempdiff.fs<- mean.fagsyl90s$meanmat[which(mean.fagsyl90s$cc=="2000-2010")]-
    mean.fagsyl90s$meanmat[which(mean.fagsyl90s$cc=="1950-1960")]
tempdiff90s.fs <- mean.fagsyl90s$meanmat[which(mean.fagsyl90s$cc=="1990-2000")]-
    mean.fagsyl90s$meanmat[which(mean.fagsyl90s$cc=="1950-1960")]
tempdiffplot.fs <- c(0, tempdiff90s.fs, tempdiff.fs)

cexhere <- 0.95
pdf(file.path("figures/pepfagsylw90s.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.5, 3), ylim=c(-6, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.fagsyl90s$cc))){
  pos.x <- tempdiffplot.fs[i]
  pos.y <- mean.fagsyl90s$matslope[i]
  sdhere <- sd.fagsyl90s$matslope[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.fagsyl90s$cc))){
  pos.x <- tempdiffplot.fs[i]
  pos.y <- mean.fagsyl90s$matslopelog[i]
  sdhere <- sd.fagsyl90s$matslopelog[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
# par(xpd=TRUE) # so I can plot legend outside
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Simple linear regression", "Using logged variables"),
   cex=1, bty="n")
dev.off()


## Now with the moving windows ...

mean.betpenswa <- aggregate(dfswa[c("matslope", "matslopelog", "meanmat")], dfswa["cc"], FUN=mean)
sd.betpenswa <- aggregate(dfswa[c("matslope", "matslopelog", "meanmat")], dfswa["cc"],  FUN=sd)

tempdiffswa <- mean.betpenswa$meanmat[which(mean.betpenswa$cc=="post")]-
    mean.betpenswa$meanmat[which(mean.betpenswa$cc=="pre")]
tempdiffplotswa <- c(0, tempdiffswa)

cexhere <- 0.95
pdf(file.path("figures/basicpepswa.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-3, 1), ylim=c(-8, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpenswa$cc))){
  pos.x <- tempdiffplotswa[i]
  pos.y <- mean.betpenswa$matslope[i]
  sdhere <- sd.betpenswa$matslope[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.betpenswa$cc))){
  pos.x <- tempdiffplotswa[i]
  pos.y <- mean.betpenswa$matslopelog[i]
  sdhere <- sd.betpenswa$matslopelog[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
# par(xpd=TRUE) # so I can plot legend outside
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Simple linear regression", "Using logged variables"),
   cex=1, bty="n")
dev.off()
