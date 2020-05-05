## Started 28 December 2019 ##
## By Lizzie ##

## Simulations of winter to spring temperatures with fstar (GDD) required increasing when chilling is low ##

# This code requires higher forcing (aka GDD, aka fstar) for leafout when chilling is low #
# It does not vary when GDD accumulates though, it always just starts accumulating on a certain day (daystostart) #

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# set.seed(113)

# Setting working directory. 
if(length(grep("ailene", getwd()))>0) { 
  setwd("~/Documents/GitHub/decsens/analyses")
} else
  setwd("~/Documents/git/projects/treegarden/decsens/analyses")


# Step 1: Set up years, days per year, temperatures, required GDD (fstar), required chill (cstar) and how much higher fstar is when cstar is not met
daysperseason <- 100
daysperinterseason <- 25
daystostart <- daysperseason+daysperinterseason # this defines the break between 'winter' and 'spring,' functionally we accumulate chill only in 1:daystostart and GDD only in daystostart:end(df)
yearz <- 20
sitez <- 45 # aka reps
degreez <- c(0, 0.5, 1, 2, 4, 6, 7) # warming -- applied evenly across whole period
sigma <- 4
fstar <- 200
cstar <- 110
fstaradjforchill <- 3 # how much more GDD to you need based on diff from cstar at end of daystostart

if(FALSE){
# Saving some intermediate progress .... setting up the seasonal temps and looking at chill and GDD calculated from them 
df <- data.frame(degwarm=numeric(), rep=numeric(), chill=numeric(), gdd=numeric())

for (i in degreez){
   for (j in 1:sitez){
       yearly_expected_temp <- rep(6, yearz)
       daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, 0 + i, sigma),
           rnorm(daysperinterseason, 2 + i , sigma), rnorm(daysperinterseason, 4 + i, sigma),
           rnorm(daysperseason, 6 + i, sigma)))
       chill <- daily_temp
       chill[(chill)<0] <- 0
       chill[(chill)>5] <- 0
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       dfadd <- data.frame(degwarm=i, rep=j, chill=mean(colSums(chill[1:daystostart,], na.rm = TRUE)),
           gdd=mean(colSums(gdd[daystostart:nrow(gdd),], na.rm = TRUE)))
       df <- rbind(df, dfadd)
    }
}

plot(daily_temp[,1]~c(1:nrow(daily_temp)), type="l")
par(mfrow=c(1,2))
plot(chill~degwarm, df)
abline(h=cstar, col="red", lwd=2)
plot(gdd~degwarm, df)
par(mfrow=c(1,1))

# Saving some more intermediate progress ... this sets up a varying fstar when chilling is below cstar
gddreq <- c()
leafout_date <- c()
for (k in 1:ncol(chill)){
    chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
    gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[daystostart:nrow(gdd),x])))
    if (chillsum[daystostart,k]>cstar) {
        gddreq[k] <- fstar
        } else {
        gddreq[k] <- fstar + (cstar-chillsum[daystostart,k])*fstaradjforchill
        }
    leafout_date[k] <- min(which(gddsum[,k] > gddreq[k]))
    }

plot(leafout_date~colMeans(daily_temp))
plot(leafout_date~gddreq)
}

## Step 2: Now I put together the seasonal temps, varying fstar (increases when chill is low) and calculate the sensitivities

df <- data.frame(degwarm=numeric(), rep=numeric(), chill=numeric(), fstar=numeric(), simplelm=numeric(),
    loglm=numeric(), perlm=numeric(),propryrschillmet=numeric(),meangddsum=numeric())#add times below cstar- times chilling not met

# yearlytemp <- "postwinter"
yearlytemp <- "alltemps"
#plot the simulated data and simple linear models from one site
site2plot = 11#arbirarily pick a site to plot! 
figname<-paste("figures/simsiteplots/decsensplot_warm","_site",site2plot,".pdf",sep="_")
pdf(figname, width = 12, height = 3)
par(mfrow=c(1, length(degreez)))
for (i in degreez){
   for (j in 1:sitez){
       yearly_expected_temp <- rep(6, yearz)#what is intention with this?
       daily_temp <- sapply(yearly_expected_temp, function(x) c(rnorm(daysperseason, 0 + i, sigma),
           rnorm(daysperinterseason, 2 + i , sigma), rnorm(daysperinterseason, 4 + i, sigma),
           rnorm(daysperseason, 6 + i, sigma)))#yearly_expected_temp not referenced at all in generating daily temps- is this intentional?
       chill <- daily_temp
       chill[(chill)<0] <- 0
       chill[(chill)>5] <- 0
       gdd <- daily_temp
       gdd[(gdd)<0] <- 0
       gddreq <- c()
       leafout_date <- c()
       for (k in 1:ncol(chill)){
           chillsum <- sapply(1:ncol(chill), function(x) (cumsum(chill[,x])))
           gddsum <- sapply(1:ncol(gdd), function(x) (cumsum(gdd[daystostart:nrow(gdd),x])))
           if (chillsum[daystostart,k]>cstar) {
           gddreq[k] <- fstar
           } else {
           gddreq[k] <- fstar + (cstar-chillsum[daystostart,k])*fstaradjforchill
           }
           leafout_date[k] <- min(which(gddsum[,k] > gddreq[k]))
           meanchill <- mean(chillsum[daystostart,])#why taking mean here? mean across 30 years?
           meanfstar <- mean(gddreq)
           chillmet<-length(which(chillsum[daystostart,]>cstar))/yearz
           meangddsum<- mean(gddsum)
           }
           yearly_tempall <- colMeans(daily_temp)
           yearly_temppostwinter <- colMeans(daily_temp[daystostart:nrow(daily_temp),])
           if(yearlytemp=="postwinter") {
               yearly_temp <- yearly_temppostwinter
               } else {
               yearly_temp <- yearly_tempall
               }
           per_leafout_date <- leafout_date/mean(leafout_date)
           per_yearly_temp <- yearly_temp/mean(yearly_temp)
           #plot(yearly_temp, leafout_date, pch=20)
           dfadd <- data.frame(degwarm=i, rep=j, chill=meanchill, fstar=meanfstar,     
               simplelm=coef(lm(leafout_date~yearly_temp))[2],
               loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
               perlm=coef(lm(per_leafout_date~per_yearly_temp))[2], propryrschillmet = chillmet,meangddsum= meangddsum)
                
                if(j==site2plot){
                  plot(yearly_temp, leafout_date,pch=16,col="darkgreen", bty="l",xlab="Temperature",ylab="Leafout", main = paste(i,"warming_site",j, sep=""))
                  m<-lm(leafout_date~yearly_temp)
                  if(summary(m)$coef[2,4]<=0.1){abline(m,lwd=2,col="darkgreen")}
                }
           df <- rbind(df, dfadd)
           
   }
   }
dev.off()
plot(simplelm~degwarm, data=df, pch=16, ylab="Sensitivity (days/C or log(days)/log(C)", xlab="Degree warming")
points(loglm~degwarm, data=df, col="dodgerblue")
points(perlm~degwarm, data=df, col="firebrick")

plot(propryrschillmet~degwarm, data=df, pch=16, ylab="Number of years (out of 30) when chilling is met", xlab="Degree warming")
plot(fstar~degwarm, data=df, pch=16, ylab="GDD", xlab="Degree warming")
plot(meangddsum~degwarm, data=df, pch=16, ylab="GDD", xlab="Degree warming")

# add cstar to df so that we can show when chilling is below cstar
#df$cstar<-cstar
#lowchill<-df$chill

##############
## Plotting ##
##############

mean.sims <- aggregate(df[c("simplelm", "loglm", "perlm","propryrschillmet", "fstar","meangddsum")], df["degwarm"], FUN=mean)
sd.sims <- aggregate(df[c("simplelm", "loglm", "perlm","propryrschillmet", "fstar","meangddsum")], df["degwarm"], FUN=sd)

cexhere <- 0.95
pdf(file.path("figures/shiftingcuessims_2panels.pdf"), width = 6, height = 8)
#par(xpd=FALSE)
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
   cex=1, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(0, 1),mgp=c(1.5,.5,0), tck=-.01,xaxs="i",yaxs = "i",
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

plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim= c(200,300),yaxt="n", ylab="",xaxt="n", xlab="", bty="u",mgp=c(1.5,.5,0), tck=-.01)
axis(side = 4,mgp=c(1.5,.5,0), tck=-.01)
mtext(expression(paste("Thermal sum required for leafout (", degree, "C)"), sep=""),side=4, adj=.5, line=2)

for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$fstar[i]
  sdhere <- sd.sims$fstar[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkred")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkred")
}
legend("topright", pch=c(19, 19), col=c("gray", "darkred"), legend=c("Chilling", "Forcing"),
       cex=1, bty="n")

dev.off()

head(df)
