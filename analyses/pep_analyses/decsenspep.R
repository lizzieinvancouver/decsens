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


##############
## Plotting ##
##############

mean.betpen <- aggregate(df[c("matslope", "matslopelog", "meanmat")], df["cc"], FUN=mean)
sd.betpen <- aggregate(df[c("matslope", "matslopelog", "meanmat")], df["cc"],  FUN=sd)

tempdiff <- mean.betpen$meanmat[which(mean.betpen$cc=="2000-2010")]-
    mean.betpen$meanmat[which(mean.betpen$cc=="1950-1960")]
tempdiffplot <- c(0, tempdiff)

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
# par(xpd=TRUE) # so I can plot legend outside
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Simple linear regression", "Using logged variables"),
   cex=1, bty="n")
dev.off()
