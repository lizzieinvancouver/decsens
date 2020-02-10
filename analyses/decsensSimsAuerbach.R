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
setwd("~/Documents/git/projects/treegarden/decsens/analyses")
plot(log(yearly_temp_trunc), log(leafout_date), pch=20, col = "dodgerblue") 
pdf(file.path("figures/simslogging.pdf"), width = 9, height = 5)
par(mfrow=c(2,3))
plot(yearly_temp, leafout_date, pch=20, main="Spring temp this row")
plot(yearly_temp, log(leafout_date), pch=20)
plot(log(yearly_temp), log(leafout_date), pch=20)
plot(yearly_temp_trunc, leafout_date, pch=20, main="Temp to leafout this row")
plot(yearly_temp_trunc, log(leafout_date), pch=20)
plot(log(yearly_temp_trunc), log(leafout_date), pch=20)
dev.off()
