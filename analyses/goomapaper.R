## Started 3 February 2021 ##
## By Lizzie ##

## See also projects/misc/pep725/pep725check_decsens.R ##

## housekeeping 
rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/projects/treegarden/decsens/analyses") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("hereliesboomboom")

##
goo <- read.csv("..//_DOthis_FLSrefs/Data.csv")
goo$X <- NULL

## plotting to find sites
whichsites1 <- subset(goo, select=c("Species", "PEP_ID", "lat", "lon"))
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


## looking at LO versus FLO over time

table(goo$Species)
table(goo$BBCH_leaf, goo$Species)
table(goo$BBCH_flower, goo$Species)

zz <- data.frame(latbi=character(), pepid=numeric(), lat=numeric(), lon=numeric(), lmslopebb=numeric(), logslopebb=numeric(), lmslopeflo=numeric(), logslopeflo=numeric(), nyears=numeric()) # bbchleaf=numeric(), bbchflo=numeric(),

subber <- subset(goo, Species=="Aesculus hippocastanum")
subberearly <- subset(subber, Year<1980)
subberearlynoNA <- subset(subberearly, is.na(Tave_yearly)==FALSE)

sitez <- unique(subber$PEP_ID)
for (site in sitez){
    subby <- subberearlynoNA[which(subberearlynoNA$PEP_ID==site),]
    if(nrow(subby)>0) {
    lmslopebb <- coef(lm(DOY_leaf~Tave_yearly, data=subby))[2]
    lmslopeflo <- coef(lm(DOY_flower~Tave_yearly, data=subby))[2]
    logslopebb <- coef(lm(log(DOY_leaf)~log(Tave_yearly), data=subby))[2]
    logslopeflo <- coef(lm(log(DOY_flower)~log(Tave_yearly), data=subby))[2]
    zzadd <- data.frame(latbi=subby$Species[1], pepid=site,
        lat=subby$lat[1], lon=subby$lon[1],
        lmslopebb=lmslopebb, logslopebb=logslopebb,
        lmslopeflo=lmslopeflo, logslopeflo=logslopeflo,
        nyears=nrow(subby))
    zz <- rbind(zz, zzadd)
    }
    }

yy <- zz

zz <- data.frame(latbi=character(), pepid=numeric(), lat=numeric(), lon=numeric(), lmslopebb=numeric(), logslopebb=numeric(), lmslopeflo=numeric(), logslopeflo=numeric(), nyears=numeric()) # bbchleaf=numeric(), bbchflo=numeric(),

subber <- subset(goo, Species=="Aesculus hippocastanum")
subberlate <- subset(subber, Year>1979)
subberlatenoNA <- subset(subberlate, is.na(Tave_yearly)==FALSE)

sitez <- unique(subber$PEP_ID)
for (site in sitez){
    subby <- subberlatenoNA[which(subberlatenoNA$PEP_ID==site),]
    if(nrow(subby)>0) {
    lmslopebb <- coef(lm(DOY_leaf~Tave_yearly, data=subby))[2]
    lmslopeflo <- coef(lm(DOY_flower~Tave_yearly, data=subby))[2]
    logslopebb <- coef(lm(log(DOY_leaf)~log(Tave_yearly), data=subby))[2]
    logslopeflo <- coef(lm(log(DOY_flower)~log(Tave_yearly), data=subby))[2]
    zzadd <- data.frame(latbi=subby$Species[1], pepid=site,
        lat=subby$lat[1], lon=subby$lon[1],
        lmslopebb=lmslopebb, logslopebb=logslopebb,
        lmslopeflo=lmslopeflo, logslopeflo=logslopeflo,
        nyears=nrow(subby))
    zz <- rbind(zz, zzadd)
    }
    }

yy$when <- "before1980"
zz$when <- "after1979"

subset(yy, pepid=="1002")
subset(zz, pepid=="1002")


par(mfrow=c(2,4))
xlimhere <- c(-20,5)
xlimherelog <- c(-1,1)
breakz <- 1000
hist(yy$lmslopebb, breaks=breakz, main="Pre CC, BB", xlim=xlimhere)
hist(zz$lmslopebb, breaks=breakz, main="Post CC, BB", xlim=xlimhere)
hist(yy$logslopebb, breaks=breakz, main="Pre CC, BB log", xlim=xlimherelog)
hist(zz$logslopebb, breaks=breakz,  main="Post CC, BB log", xlim=xlimherelog)
hist(yy$lmslopeflo, breaks=breakz, main="Pre CC, Flo", xlim=xlimhere)
hist(zz$lmslopeflo, breaks=breakz, main="Post CC, Flo", xlim=xlimhere)
hist(yy$logslopeflo, breaks=breakz, main="Pre CC, Flo log", xlim=xlimherelog)
hist(zz$logslopeflo, breaks=breakz,  main="Post CC, Flo log", xlim=xlimherelog)

mean(yy$lmslopebb, na.rm=TRUE)
mean(zz$lmslopebb, na.rm=TRUE)
mean(yy$logslopebb, na.rm=TRUE)
mean(zz$logslopebb, na.rm=TRUE)

mean(yy$lmslopeflo, na.rm=TRUE)
mean(zz$lmslopeflo, na.rm=TRUE)
mean(yy$logslopeflo, na.rm=TRUE)
mean(zz$logslopeflo, na.rm=TRUE)
