## Started 19 April 2020 ##
## By Lizzie ##

## Look at the Betpen data for signal of non-linearities ... 

# housekeeping
rm(list=ls()) 
options(stringsAsFactors=FALSE)

setwd("~/Documents/git/projects/treegarden/decsens/analyses")

# get the data from ospree repo
if(FALSE){
osp <- read.csv("..//..//budreview/ospree/analyses/output/ospree_clean_withchill_BB.csv", header = TRUE)
bpall <- subset(osp, genus=="Betula" & species=="pendula")
source("..//..//budreview/ospree/analyses/source/commoncols.R")
bp <- bpall[,which(colnames(bpall) %in% c(common.cols.wchill, "forcetemp_night", "photoperiod_night"))]

bp$forceday <- as.numeric(bp$forcetemp)
bp$forcenight <- as.numeric(bp$forcetemp_night)
bp$photonight <- as.numeric(bp$photoperiod_night)

bp$photo <- as.numeric(bp$photoperiod_day)
bp$force <- bp$forceday
# Adjust forcing temperature to incorporate day/night temp differentials, if present
bp$force[is.na(bp$forcenight)==FALSE & is.na(bp$photo)==FALSE &
    is.na(bp$photonight)==FALSE] <-
    (bp$forceday[is.na(bp$forcenight)==FALSE & is.na(bp$photo)==FALSE &
    is.na(bp$photonight)==FALSE]*
    bp$photo[is.na(bp$forcenight)==FALSE & is.na(bp$photo)==FALSE &
    is.na(bp$photonight)==FALSE] +
    bp$forcenight[is.na(bp$forcenight)==FALSE & is.na(bp$photo)==FALSE &
    is.na(bp$photonight)==FALSE]*
    bp$photonight[is.na(bp$forcenight)==FALSE & is.na(bp$photo)==FALSE &
    is.na(bp$photonight)==FALSE])/24

write.csv(bp, "cherry/betpen_ospree.csv", row.names=FALSE)
}


bp <- read.csv("cherry/betpen_ospree.csv")
bp$resp <- as.numeric(bp$response.time)
bp$gddreq <- bp$force * bp$resp

# plot required thermal sum against one metric (Utah) of chilling
library(ggplot2)
# plot, removing junttilla which is a dormancy release study
ggplot(subset(bp, Total_Utah_Model>0), aes(Total_Utah_Model, gddreq, col=as.factor(photo))) +
    geom_point() + facet_wrap(datasetID ~ .)
