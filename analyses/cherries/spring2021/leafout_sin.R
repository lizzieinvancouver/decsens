library("tidyverse")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

################################
## Accumulation date "a" known #
################################

##Simulated Data

#set unknown parameters
alpha_0 <- 60      # average annual temperature
alpha_1 <- 30      # temperature amplitude (max/min temp is 60 +/- 30)
alpha_2 <- 90      # phase in days (90 means min temp occurs on day 0)
sigma   <- 20      # temperature sd (any day is off the average by 20 degrees)
s       <- 3000    # leafout happens when cumulative temperatures reach s

#set known parameters
a       <- 9       # day plant starts accumulating temperature
b       <- 100     # last day of temperature window 
                   ## (i.e. average temperature taken from day 0 to b)
n_sim   <- 100     # number of experiments

#function to generate expected daily temperature
daily_temp <- function(n) alpha_0 + alpha_1 * sin(2 * pi / 365 * (n - alpha_2))

#create data
data <- data.frame(leaf_date = numeric(0),
                   mean_temp = numeric(0))

for(sim in 1:n_sim) {
  temp      <- sapply(1:b, daily_temp) + rnorm(b, 0, sigma)
  leaf_date <- which.min(cumsum(temp[(a+1):b]) < s) + a
  mean_temp <- mean(temp)
  data      <- rbind(data, data.frame(leaf_date, mean_temp))
}

##Model Simulated Data
stan_data <- list(N       = n_sim,
                  M       = data$mean_temp,
                  n_s     = data$leaf_date,
                  a       = a,
                  b       = b)

init <- list(
  list(alpha_0 = 160, alpha_1 = 40, alpha_2 = 100, s = 1000, sigma = 1000),
  list(alpha_0 = 80,  alpha_1 = 20, alpha_2 = 100, s = 1000, sigma = 1000),
  list(alpha_0 = 40,  alpha_1 = 10, alpha_2 = 100, s = 1000, sigma = 1000),
  list(alpha_0 = 20,  alpha_1 = 5,  alpha_2 = 100, s = 1000, sigma = 1000)
)

model <- stan_model("leafout_sin_a.stan")
mle <- optimizing(model, data = stan_data, init = init[[1]], seed = 1)
fit <- sampling(model, data = stan_data, init = init, seed = 1, chains = 4)

####################################
## Accumulation date "a" not known #
####################################

##Simulated Data

#set unknown parameters
alpha_0 <- 60      # average annual temperature
alpha_1 <- 30      # temperature amplitude (max/min temp is 60 +/- 30)
alpha_2 <- 90      # phase in days (90 means min temp occurs on day 0)
sigma   <- 20      # temperature sd (any day is off the average by 20 degrees)
s       <- 3000    # leafout happens when cumulative temperatures reach s
p       <- .1      # probability plant starts collecting temperature given 
                   ## it hasn't already started collecting temperature

#set known parameters
b       <- 100     # last day of temperature window 
                   ## (i.e. average temperature taken from day 0 to b)
n_sim   <- 100     # number of experiments

#function to generate expected daily temperature
daily_temp <- function(n) alpha_0 + alpha_1 * sin(2 * pi / 365 * (n - alpha_2))

#create data
data <- data.frame(leaf_date = numeric(0),
                   mean_temp = numeric(0))

for(sim in 1:n_sim) {
  a         <- 1 + rgeom(1, p) # day plant starts accumulating temperature
                               ## i.e. plant collects from day a to leaf_date
  temp      <- sapply(1:b, daily_temp) + rnorm(b, 0, sigma)
  leaf_date <- which.min(cumsum(temp[a:b]) < s) + a
  mean_temp <- mean(temp)
  data      <- rbind(data, data.frame(leaf_date, mean_temp))
  }

##Model Simulated Data
stan_data <- list(N       = n_sim,
                  M       = data$mean_temp,
                  n_s     = data$leaf_date,
                  b       = b,
                  alpha_1   = alpha_1)

init <- list(
  list(alpha_0 = 160, alpha_1 = 40, alpha_2 = 100, s = 1000, sigma = 1000, p = .1),
  list(alpha_0 = 80,  alpha_1 = 20, alpha_2 = 100, s = 1000, sigma = 1000, p = .1),
  list(alpha_0 = 40,  alpha_1 = 10, alpha_2 = 100, s = 1000, sigma = 1000, p = .1),
  list(alpha_0 = 20,  alpha_1 = 5,  alpha_2 = 100, s = 1000, sigma = 1000, p = .1)
)

model <- stan_model("leafout_sin_p.stan")
#mle <- optimizing(model, data = stan_data, init = init[[1]], seed = 1) 
## LBFGS doesn't seem to work with this model
fit <- sampling(model, data = stan_data, init = init, seed = 1, chains = 4)