## Started 3 February 2021 ##
## By Lizzie ##

## Data (PEP725) from Ma et al. 2020 (GCB) showing first events are more sensitive (duh)
# "We found that regardless whether flowering or leaf-out occurred first, the first event advanced more than the second during 1950–2013, resulting in a prolonged time interval between the two events"
# So, we want to see how similar the slopes become if you log the events and temperature ...

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

## libraries
library(ggplot2)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/projects/treegarden/decsens/analyses") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("hereliesboomboom")

##
madat <- read.csv("ma2020/Data.csv")
madat$X <- NULL

## plotting to get a sense sites
whichsites1 <- subset(madat, select=c("Species", "PEP_ID", "lat", "lon"))
whichsites <- whichsites1[!duplicated(whichsites1),]

library("maptools")
data(wrld_simpl)
EuropeList <- c('Germany', 'France', 'Sweden', 'Finland', 'Norway', 'Austria', 'Switzerland')
my_map <- wrld_simpl[wrld_simpl$NAME %in% EuropeList, ]
plot(my_map)

whichsitesAH <- subset(whichsites, Species=="Aesculus hippocastanum")
whichsitesAG <- subset(whichsites, Species=="Alnus glutinosa")
whichsitesFE <- subset(whichsites, Species=="Fraxinus excelsior")

points(whichsitesAH$lon, whichsitesAH$lat, pch=16, col="red")
points(whichsitesAG$lon, whichsitesAG$lat, pch=16, col="darkgreen")
points(whichsitesFE$lon, whichsitesFE$lat, pch=16, col="dodgerblue") 

table(madat$Species)
table(madat$BBCH_leaf, madat$Species)
table(madat$BBCH_flower, madat$Species)

zz <- data.frame(latbi=character(),
                 pepid=numeric(),
                 lat=numeric(),
                 lon=numeric(),
                 lmslopebb=numeric(),
                 logslopebb=numeric(),
                 lmslopeflo=numeric(),
                 logslopeflo=numeric(),
                 nyears=numeric(),
                 avgTemp=numeric()) 


madatnoNA <- subset(madat, is.na(Tave_yearly)==FALSE)
# T_ave site is much longer than siteID ...?
length(unique(madat$Tave_site))
length(unique(madat$PEP_ID))

for (latbihere in unique(madat$Species)){
    subber <- subset(madatnoNA, Species==latbihere)
    for (sitehere in unique(subber$PEP_ID)){
        subby <- subber[which(subber$PEP_ID==sitehere),]
            if(nrow(subby)>0) {
                lmslopebb <- coef(lm(DOY_leaf~Tave_yearly, data=subby))[2]
                lmslopeflo <- coef(lm(DOY_flower~Tave_yearly, data=subby))[2]
                logslopebb <- coef(lm(log(DOY_leaf)~log(Tave_yearly), data=subby))[2]
                logslopeflo <- coef(lm(log(DOY_flower)~log(Tave_yearly), data=subby))[2]
                zzadd <- data.frame(latbi=latbihere,
                            pepid=sitehere,
                            lat=subby$lat[1],
                            lon=subby$lon[1],
                            lmslopebb=lmslopebb,
                            logslopebb=logslopebb,
                            lmslopeflo=lmslopeflo,
                            logslopeflo=logslopeflo,
                            nyears=nrow(subby),
                            avgTemp=mean(subby$Tave_site, na.rm=TRUE))
                zz <- rbind(zz, zzadd)
                }
        }
    }


aggregate(zz[c("lmslopebb", "logslopebb", "lmslopeflo", "logslopeflo")], zz[c("latbi")], FUN=mean)
# Pretty similar results with lm and log ...
# but it's identical climate so I think it would have to be, right?

bbzz <- subset(zz, select=c("latbi", "pepid", "lat", "lon", "lmslopebb", "logslopebb"))
names(bbzz) <- c("latbi", "pepid", "lat", "lon", "lmslope", "logslope")
bbzz$event <- "budburst"
flozz <- subset(zz, select=c("latbi", "pepid", "lat", "lon", "lmslopeflo", "logslopeflo"))
names(flozz) <- c("latbi", "pepid", "lat", "lon", "lmslope", "logslope")
flozz$event <- "flowering"

zzlong <- rbind(bbzz, flozz)

ggplot(zzlong, aes(x=event, y=lmslope)) +
    geom_boxplot() +
    facet_wrap(.~latbi, scales="free") +
    xlab("Event") +
    ylab("lm slope") +
    theme_minimal()

ggplot(zzlong, aes(x=event, y=logslope)) +
    geom_boxplot() +
    facet_wrap(.~latbi, scales="free") +
    xlab("Event") +
    ylab("log slope") +
    theme_minimal()

# By latitude ...
zz45 <- subset(zz, lat>45)

ggplot(zz45, aes(x=lat, y=lmslopebb)) +
    geom_point() +
    geom_smooth(method = "lm", linetype = 2, lwd=0.5, color="dodgerblue", se = FALSE) +
    facet_wrap(.~latbi, scales="free") +
    xlab("Event") +
    ylab("lm slope") +
    theme_minimal()

quartz()
ggplot(zz45, aes(x=lat, y=logslopebb)) +
    geom_point() +
    geom_smooth(method = "lm", linetype = 2, lwd=0.5, color="dodgerblue", se = FALSE) +
    facet_wrap(.~latbi, scales="free") +
    xlab("Event") +
    ylab("log slope") +
    theme_minimal()

