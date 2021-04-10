library("tidyverse")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

##Simulated Data

#set unknown parameters
alpha_0 <- 32      # temperature at day 0
alpha_1 <- .4      # temperature increases an average of .4 degrees F each day
sigma   <- 20      # temperature sd (any day is off the average by 20 degrees)
s       <- 3000    # leafout happens when cumulative temperatures reach s
p       <- .1      # probability plant starts collecting temperature given 
                   ## it hasn't already started collecting temperature

#set known parameters
b       <- 100     # last day temperature window (average temperature taken from 
                   ## day 0 to b)
n_sim   <- 100     # number of experiments

#create data
data <- data.frame(leaf_date = numeric(0),
                   mean_temp = numeric(0))

for(sim in 1:n_sim) {
  a         <- 1 + rgeom(1, p) # the day the plant starts collecting temperature
  temp      <- alpha_0 + alpha_1 * (1:b) + rnorm(b, 0, sigma)
  leaf_date <- which.min(cumsum(temp[a:b]) < s) + a
  mean_temp <- mean(temp)
  data      <- rbind(data, data.frame(leaf_date, mean_temp))
  }

##Model Simulated Data
stan_data <- list(N       = n_sim,
                  M       = data$mean_temp,
                  n_s     = data$leaf_date,
                  b       = b)

init <- function () 
  list(alpha_0 = 1,
       alpha_1 = 1, 
       s = 1, 
       sigma = 100,
       p = .01)

model <- stan_model("leafout.stan")
fit <- sampling(model, data=stan_data, init = init, seed = 1, chains = 4,
                control = list(adapt_delta = 0.99))
