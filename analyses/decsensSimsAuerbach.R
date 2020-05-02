## Started 18 December 2019 ##
## Code by Jonathan Auerbach (and Andrew Gelman some) ##
## All comments by Lizzie ##

## From Jonathan's email (5 April 2020) ##
# "Even if the linear relationship were true, I think y ~ x is the wrong regression. The researchers stop accumulating the average x based on y. So really they should run the inverse regression lm(x ~ y) and calibrate. This won't be much of a problem in sample, but if researchers are extrapolating y for higher x, it will lead to substantial bias."

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


