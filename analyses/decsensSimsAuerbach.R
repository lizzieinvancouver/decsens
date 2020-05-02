## Started 18 December 2019 ##
## Code by Jonathan Auerbach (and Andrew Gelman some) ##
## All comments by Lizzie ##

## Quick simulation of the declining sensitivities problem ##

# First, build up some simulated climate and leafout (after 200 GDD) data
yearly_expected_temp <- c(rep(0,20), rep(5,20), rep(10,20)) # create holder for three sets of temps, 20 per set (note, these are bigger than we have done -- we have done 1 to 7 degrees, this does 0, 5 and 10 degrees)

# Or, do 4 sets ... (Lizzie added this next line)
yearly_expected_temp <- c(rep(0,20), rep(5,20), rep(10,20), rep(20, 20)) 
daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(30, 10 + x, 4)) # now make some simple daily temps
# daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(30, 10 + x, (4+x))) # Lizzie added
leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > 200))) # set leafout date as whenever 200 GDD is reached
yearly_temp <- colMeans(daily_temp) # estimate the mean temp of each simulated dataset

yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x])) # estimate the mean temp of each simulated dataset only until leafout 

plot(yearly_temp, leafout_date, pch=20)
points(yearly_temp_trunc, leafout_date, pch=20, col = "red")

plot(yearly_temp_trunc, leafout_date, pch=20, col = "red")

# Lizzie added all of the below
cexhere <- 0.5
setwd("~/Documents/git/projects/treegarden/decsens/analyses")
plot(log(yearly_temp_trunc), log(leafout_date), pch=20, col = "dodgerblue") 
pdf(file.path("figures/simslogging.pdf"), width = 9, height = 5)
par(mfrow=c(2,3))
plot(yearly_temp_trunc, leafout_date, pch=20, xlab="Simulated spring temperature to leafout",
     ylab="Leafout date", main="", cex=cexhere)
plot(yearly_temp_trunc, log(leafout_date), pch=20, xlab="Simulated spring temperature to leafout",
     ylab="log(Leafout date)", main="", cex=cexhere)
plot(log(yearly_temp_trunc), log(leafout_date), pch=20, xlab="log(Simulated spring temperature to leafout)",
     ylab="log(Leafout date)", main="", cex=cexhere)
plot(yearly_temp, leafout_date, pch=20, xlab="Simulated spring temperature",
    ylab="Leafout date", main="", cex=cexhere)
plot(yearly_temp, log(leafout_date), pch=20, xlab="Simulated spring temperature",
    ylab="log(Leafout date)", main="", cex=cexhere)
plot(log(yearly_temp), log(leafout_date), pch=20, xlab="log(Simulated spring temperature)",
     ylab="log(Leafout date)", main="", cex=cexhere)
dev.off()


## From Jonathan's email (5 April 2020) ##

#part 1
data <- data.frame(leaf_date = numeric(0),
                   cum_temp  = numeric(0),
                   mean_temp = numeric(0),
                   threshold = numeric(0),
                   delta     = numeric(0))

threshold <- 1000
for(delta in c(5, 10, 15, 20)) {
  for(sim in 1:1000) {
    temp <- delta * (1:100) + rnorm(100, 0, 50)
    leaf_date <- which.min(cumsum(temp) < threshold)
    cum_temp <- sum(temp[1:leaf_date])
    mean_temp <- mean(temp[1:leaf_date])
    data <- rbind(data, data.frame(leaf_date, cum_temp, mean_temp, threshold, delta))
  }
}

sapply(unique(data$delta), function(i) sd(data$leaf_date[data$delta == i]))
MASS::boxcox(lm(leaf_date ~ cum_temp, data))

#part 2
library(ggplot2)
yxplot<- ggplot() +
  aes(mean_temp, leaf_date) + # leaf ~ temp
  geom_point(data = data) +
  theme_bw() +
  xlab("Mean temperature") + ylab("Leafout date") +
  geom_smooth(method = "lm", color = "blue", fullrange = TRUE,
              data = data[data$mean_temp < 100, ]) +
  geom_smooth(method = "glm", formula = y~log(x), color = "red",
              fullrange=TRUE,
              method.args = list(family = gaussian(link = 'log')),
              data = data[data$mean_temp < 100, ]) +
  coord_cartesian(ylim = c(5, 30), xlim = c(30, 160)) # coord_cartesian

xyplot <- ggplot() +
  aes(leaf_date, mean_temp) + # temp ~ leaf
  geom_point(data = data) +
  theme_bw() +
  ylab("Mean temperature") + xlab("Leafout date") +
  geom_smooth(method = "lm", color = "blue", fullrange = TRUE,
              data = data[data$mean_temp < 100, ]) +
  geom_smooth(method = "glm", formula = y~log(x), color = "red",
              fullrange=TRUE,
              method.args = list(family = gaussian(link = 'log')),
              data = data[data$mean_temp < 100, ]) +
  coord_flip(xlim = c(5, 30), ylim = c(30, 160)) # coord_flip



# Setting working directory. 
if(length(grep("ailene", getwd()))>0) { 
  setwd("~/Documents/GitHub/decsens/analyses")
} else
  setwd("~/Documents/git/projects/treegarden/decsens/analyses")

library(gridExtra)
plotsave <- grid.arrange(yxplot, xyplot, nrow = 1)

ggsave("figures/compareyxxy.pdf", plot=plotsave, width = 7, height = 4)


