## Started 17 January 2021 ##
## By Lizzie ##
## Happy Indoguration! ##

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/projects/treegarden/decsens/analyses") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("hereliesboomboom")


library(ggplot2)

# Set up time periods and noise
dayz <- 366
yearz <- 5
sigma_y <- 5

# Set up one sine curve
amplitude =  20
center = 60
width =  dayz/2
offyoset =  10 # controls y axis

x <- (1:(dayz*yearz))
y = amplitude * sin(pi * (x - center) / width) + offyoset
ynoisy <- y + rnorm(length(y), 0, sigma_y)
plot(ynoisy~x, type="l")

# yalt = amplitude * sin(2 * pi * (x - center) / width) + offyoset # Geoff suggested this but it makes no sense to me
# plot((yalt + rnorm(length(y), 0, sigma_y))~x, type="l", main="added the 2 as suggested")

df <- data.frame(doy=rep(1:dayz, yearz), clim=ynoisy, year=rep(1:yearz, each=dayz))
ggplot(df, aes(x=doy, y=clim)) +
    geom_line(color="dodgerblue") +
    facet_wrap(.~as.factor(year)) +
    xlab("Day") +
    ylab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()

## Playing around with offset
offyosets <- c(0, 5, 10, 15, 20, 30)
par(mfrow=c(3, 2))
for (offset in offyosets){
    y = amplitude * sin(pi * (x - center) / width) + offset
    ynoisy <- y + rnorm(length(y), 0, sigma_y)
    plot(ynoisy~x, type="l", main=paste("offset", offset))
    }

##
## Okay, climate done! Onward to leafout
yearz <- 50
threshold <- 200
biozero <- -2 # above what temperature do plants accumulate

x <- (1:(dayz*yearz))
y = amplitude * sin(pi * (x - center) / width) + offyoset
ynoisy <- y + rnorm(length(y), 0, sigma_y)
plot(ynoisy~x, type="l")
df <- data.frame(doy=rep(1:dayz, yearz), clim=ynoisy, year=rep(1:yearz, each=dayz))

# Just get leafout across years ... (an example really)
dfleafy <- data.frame(year=numeric(), leafout_date=numeric(), annual_temp=numeric(), threemon_temp=numeric())
for (i in 1:yearz){
    subby <- subset(df, year==i)
    leafout_date <- min(which(cumsum(subby$clim) > threshold))
    annual_temp <- mean(subby$clim)
    threemon_temp <- mean(subby$clim[30:120])
    dfadd <- data.frame(year=i, leafout_date=leafout_date, annual_temp=annual_temp, threemon_temp=threemon_temp)
    dfleafy <- rbind(dfleafy, dfadd) 
    }

## Some sims across climate space...
# Vary just the offset for now 
offyosets <- c(0, 5, 10, 15, 20, 30)

# The below does offset x year inefficiently 
dfleafyslow <- data.frame(year=numeric(), leafout_date=numeric(), annual_temp=numeric(), threemon_temp=numeric(),
    offset=numeric())
for (offset in offyosets){
    y = amplitude * sin(pi * (x - center) / width) + offset
    ynoisy <- y + rnorm(length(y), 0, sigma_y)
    df <- data.frame(doy=rep(1:dayz, yearz), clim=ynoisy, year=rep(1:yearz, each=dayz))
    for (i in 1:yearz){
        subby <- subset(df, year==i)
        dfadd <- data.frame(year=i, leafout_date=min(which(cumsum(subby$clim) > threshold)), 
            annual_temp=mean(subby$clim), threemon_temp=mean(subby$clim[30:120]), offset=offset)
        dfleafyslow <- rbind(dfleafyslow, dfadd) 
    }
}

ggplot(dfleafyslow, aes(x=annual_temp, y=leafout_date)) +
    geom_point(color="dodgerblue") +
    geom_smooth(method = "lm", linetype = 2, lwd=0.5, color="darkgray", se = FALSE) +
    facet_wrap(.~as.factor(offset)) +
    ylab("Leafout day") +
    xlab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()

# The below does offset x year more efficiently
dfleafy <- data.frame(leafout_date = numeric(0),
                   cum_temp  = numeric(0),
                   mean_temp = numeric(0),
                   threemon_temp = numeric(0),
                   threshold = numeric(0),
                   offset = numeric(0))

for (offset in offyosets){
    for (i in 1:yearz){
        y = amplitude * sin(pi * (x - center) / width) + offset
        ybio <- y + rnorm(length(y), 0, sigma_y)
        ybio[which(ybio<biozero)] <- biozero
        leafout_date <- which.min(cumsum(ybio) < threshold)
        cum_temp <- sum(ybio[1:leafout_date])
        mean_temp <- mean(ybio[1:leafout_date])
        threemon_temp <- mean(ybio[30:120])
        dfleafy <- rbind(dfleafy, data.frame(leafout_date, cum_temp, mean_temp, threshold, offset))
        }
}

sapply(unique(dfleafy$offset), function(i) coef(lm(dfleafy$leafout_date[dfleafy$offset == i]~
    dfleafy$mean_temp[dfleafy$offset == i])))

ggplot(dfleafy, aes(x=mean_temp, y=leafout_date, group=offset, color=offset)) +
    geom_point() +
    geom_smooth(method = "lm", linetype = 2, lwd=0.5, color="darkgray", se = FALSE) +
    # facet_wrap(.~as.factor(offset)) +
    ylab("Leafout day") +
    xlab(expression(paste("Daily temperature (", degree, "C)"), sep="")) +
    theme_minimal()

if(FALSE){
## To do  ...
# Add in fall date dependent on just photoperiod with
library(geosphere) # for daylengths
photoper <- daylength(45, 1:dayz) # latitude=45
pstar <- 14 # threshold to leafout
pstarday <- max(which(photoper > pstar))
}
