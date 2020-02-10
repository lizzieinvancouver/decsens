## Started 18 January 2020 ##
## By Lizzie ##

## PEP data ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Setting working directory. 
setwd("~/Documents/git/projects/treegarden/decsens/analyses")

# Get data ...
df <- read.csv("pep_analyses/output/bpenestimates_withlog_1950-2010.csv", header=TRUE)
dfswa <- read.csv("pep_analyses/output/swaestimates_withlog.csv", header=TRUE)
fs <- read.csv("pep_analyses/output/fsylestimates_withlog_1950-2010.csv", header=TRUE)


##############
## Plotting ##
##############

mean.betpen <- aggregate(df[c("matslope", "matslopelog", "meanmat", "varmat", "varlo", "meangdd",
    "matslopeconfint11", "matslopeconfint89", "matslopelogconfint11", "matslopelogconfint89")],
    df["cc"], FUN=mean)

tempdiff1 <- mean.betpen$meanmat[which(mean.betpen$cc=="1970-1990")]-
    mean.betpen$meanmat[which(mean.betpen$cc=="1950-1970")]
tempdiff2 <- mean.betpen$meanmat[which(mean.betpen$cc=="1990-2010")]-
    mean.betpen$meanmat[which(mean.betpen$cc=="1950-1970")]

tempdiffplot <- c(0, tempdiff1, tempdiff2)


cexhere <- 0.95
cextext <- 0.5
pdf(file.path("figures/basicpep1950to2000.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslope[i]
  ciherelo <- mean.betpen$matslopeconfint11[i]
  cihereup <- mean.betpen$matslopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen$cc)[i], cex=cextext, col="darkblue")
  }
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslopelog[i]
  ciherelo <- mean.betpen$matslopelogconfint11[i]
  cihereup <- mean.betpen$matslopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=1, bty="n")
dev.off()


## For Fagus sylvatica
mean.fs <- aggregate(fs[c("matslope", "matslopelog", "meanmat", "varmat", "varlo", "meangdd",
    "matslopeconfint11", "matslopeconfint89", "matslopelogconfint11", "matslopelogconfint89")],
    fs["cc"], FUN=mean)

tempdiff1fs <- mean.fs$meanmat[which(mean.fs$cc=="1970-1990")]-
    mean.fs$meanmat[which(mean.fs$cc=="1950-1970")]
tempdiff2fs <- mean.fs$meanmat[which(mean.fs$cc=="1990-2010")]-
    mean.fs$meanmat[which(mean.fs$cc=="1950-1970")]

tempdiffplotfs <- c(0, tempdiff1fs, tempdiff2fs)

cexhere <- 0.95
pdf(file.path("figures/basicpep1950to2000fs.pdf"), width = 6, height = 4)
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="")
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.fs$cc))){
  pos.x <- tempdiffplotfs[i]
  pos.y <- mean.fs$matslope[i]
  ciherelo <- mean.fs$matslopeconfint11[i]
  cihereup <- mean.fs$matslopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen$cc)[i], cex=cextext, col="darkblue")
  }
for(i in 1:length(unique(mean.fs$cc))){
  pos.x <- tempdiffplotfs[i]
  pos.y <- mean.fs$matslopelog[i]
  ciherelo <- mean.fs$matslopelogconfint11[i]
  cihereup <- mean.fs$matslopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=1, bty="n")
dev.off()


##
## 2 panel with FagSyl -- showing raw on top panels, and just logged on bottom panels
fagsyljitter <- 0.04
cexhere <- 0.5
cexhereleg <- 0.7
pdf(file.path("figures/basicpep1950to20002spp2panel.pdf"), width = 5, height = 7)
par(xpd=FALSE)
par(mfrow=c(2,1))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")))
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslope[i]
  ciherelo <- mean.betpen$matslopeconfint11[i]
  cihereup <- mean.betpen$matslopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.fs$cc))){
  pos.x <- tempdiffplotfs[i]
  pos.y <- mean.fs$matslope[i]
  ciherelo <- mean.fs$matslopeconfint11[i]
  cihereup <- mean.fs$matslopeconfint89[i]
  lines(x=rep(pos.x+fagsyljitter, 2), y=c(ciherelo, cihereup), col="dodgerblue")
  points(pos.x+fagsyljitter, pos.y, cex=cexhere, pch=19, col="dodgerblue")
  text(pos.x + 0.03, pos.y, labels=unique(mean.fs$cc)[i], cex=cextext, col="black")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "dodgerblue"),
   legend=c("Betula pendula", "Fagus sylvatica"), cex=cexhereleg, bty="n")

# Log-log 
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-0.9, 0.1),
     ylab=expression(paste("Estimated sensitivity (log(days)/log(", degree, "C))"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")))
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslopelog[i]
  ciherelo <- mean.betpen$matslopelogconfint11[i]
  cihereup <- mean.betpen$matslopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
for(i in 1:length(unique(mean.fs$cc))){
  pos.x <- tempdiffplotfs[i]
  pos.y <- mean.fs$matslopelog[i]
  ciherelo <- mean.fs$matslopelogconfint11[i]
  cihereup <- mean.fs$matslopelogconfint89[i]
  lines(x=rep(pos.x+fagsyljitter, 2), y=c(ciherelo, cihereup), col="pink")
  points(pos.x+fagsyljitter, pos.y, cex=cexhere, pch=19, col="pink")
  text(pos.x + 0.03, pos.y, labels=unique(mean.fs$cc)[i], cex=cextext, col="black")
  }
legend("bottomright", pch=c(19, 19), col=c("salmon", "pink"), legend=c("Betula pendula", "Fagus sylvatica"),
   cex=cexhereleg, bty="n")
dev.off()
## END: Two panel with FagSyl
## 


##
## Four panel with FagSyl -- showing logged and raw on left panels, and just logged on right panels
cexhere <- 0.5
cexhereleg <- 0.7
pdf(file.path("figures/basicpep1950to20002spp4panel.pdf"), width = 9, height = 6)
par(xpd=FALSE)
par(mfrow=c(2,2))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Betpen")
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslope[i]
  ciherelo <- mean.betpen$matslopeconfint11[i]
  cihereup <- mean.betpen$matslopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen$cc)[i], cex=cextext, col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslopelog[i]
  ciherelo <- mean.betpen$matslopelogconfint11[i]
  cihereup <- mean.betpen$matslopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"),
   legend=c("Using raw x, y", "Using logged x, y"), cex=cexhereleg, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.9, 0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Betpen (logged only)")
for(i in 1:length(unique(mean.betpen$cc))){
  pos.x <- tempdiffplot[i]
  pos.y <- mean.betpen$matslopelog[i]
  ciherelo <- mean.betpen$matslopelogconfint11[i]
  cihereup <- mean.betpen$matslopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen$cc)[i], cex=cextext, col="salmon")
  }
# FagSyl 
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5), ylim=c(-10, 1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Fagsyl")
abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.fs$cc))){
  pos.x <- tempdiffplotfs[i]
  pos.y <- mean.fs$matslope[i]
  ciherelo <- mean.fs$matslopeconfint11[i]
  cihereup <- mean.fs$matslopeconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen$cc)[i], cex=cextext, col="darkblue")
  }
for(i in 1:length(unique(mean.fs$cc))){
  pos.x <- tempdiffplotfs[i]
  pos.y <- mean.fs$matslopelog[i]
  ciherelo <- mean.fs$matslopelogconfint11[i]
  cihereup <- mean.fs$matslopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=cexhereleg, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.1, 1.5),  ylim=c(-0.9, 0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C)"), sep=""),
         xlab=expression(paste("Warming (", degree, "C)")), main="Fagsyl (logged only)")
for(i in 1:length(unique(mean.fs$cc))){
  pos.x <- tempdiffplotfs[i]
  pos.y <- mean.fs$matslopelog[i]
  ciherelo <- mean.fs$matslopelogconfint11[i]
  cihereup <- mean.fs$matslopelogconfint89[i]
  lines(x=rep(pos.x, 2), y=c(ciherelo, cihereup), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  text(pos.x + 0.1, pos.y, labels=unique(mean.betpen$cc)[i], cex=cextext, col="salmon")
  }
dev.off()
## END: Four panel with FagSyl -- showing logged and raw on left panels, and just logged on right panels
## 




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
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Using raw x, y", "Using logged x, y"),
   cex=1, bty="n")
dev.off()
