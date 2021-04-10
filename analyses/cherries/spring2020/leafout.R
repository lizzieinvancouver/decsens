library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/Documents/git/projects/treegarden/decsens/analyses/cherry")

##Simulated Data

#set parameters
alpha_0 <- 32      # plants start storing on FEB 1 when average temperature exceeds freezing 
alpha_1 <- .3      # temperature increases an average of .3 degrees F each day
a <- 30            # March 1st, beginning of window for average temperature
b <- 121           # May 31st, end of window for average temperature
sigma <- 10        # temperature sd (any day is off the average by around 10 degrees)
beta <- 1500       # leafout happens when cumulative temperatures reach threshold (expectation: March 10)
n_sim <- 1000      # number of experiments

#create data
data <- data.frame(leaf_date = numeric(0),
                   mean_temp = numeric(0),
                   beta = numeric(0),
                   delta     = numeric(0))

for(sim in 1:n_sim) {
  temp <- alpha_0 + alpha_1 * (1:b) + rnorm(b, 0, sigma)
  leaf_date <- which.min(cumsum(temp) < beta)
  mean_temp <- mean(temp[a:b])
  data <- rbind(data, data.frame(leaf_date, mean_temp, beta, alpha_1))
  }


##Model with Simulated Data
#assume alpha_0 known (freezing)
stan_data <- list(N       = n_sim,
                  M       = data$mean_temp,
                  n_beta  = data$leaf_date,
                  a       = a,
                  b       = b,
                  alpha_0 = alpha_0)

model <- stan_model("leafout.stan")
fit_bayes <- sampling(model, data=stan_data, init_r = 10, seed = 1)

fit_mle <- optimizing(model, data=stan_data, hessian = TRUE, 
                      algorithm = "Newton", seed = 1)
#MLE
fit_mle$par
#Standard Errors
c(sqrt(diag(-solve(fit_mle$hessian))),
  fit_mle$par["beta"]/fit_mle$par["beta_unit_scale"] * 
    sqrt(diag(-solve(fit_mle$hessian)))["beta_unit_scale"])


##Model with Real Data

#DC Temp: https://www.weather.gov/media/lwx/climate/dcatemps.pdf
#Cherry Blossom: https://www.epa.gov/sites/production/files/2016-08/cherry-blossoms_fig-1.csv
cherry <-
  read.csv("https://www.epa.gov/sites/production/files/2016-08/cherry-blossoms_fig-1.csv",
           skip = 6)
dc_temp <- read.csv("dc_temp.csv")
dc_temp <- dc_temp[dc_temp$YEAR %in% cherry$Year, ]

n_beta <- cherry$Yoshino.peak.bloom.date # -80
M      <- rowMeans(dc_temp[, c("FEB", "MAR", "APR", "MAY")])

stan_data <- list(N       = length(n_beta),
                  M       = M,
                  n_beta  = n_beta,
                  a       = a,
                  b       = b,
                  alpha_0 = alpha_0)

fit_bayes <- sampling(model, data=stan_data, init_r = 10, seed = 1)
fit_bayes

if(FALSE){
## check with equinox correction (does little)
library(timeDate)
cherry$verneq <- as.numeric(format(JPVernalEquinox(cherry$Year), "%j"))
cherry$Yoshino.peak.bloom.date.adj <- cherry$Yoshino.peak.bloom.date-cherry$verneq

n_beta <- cherry$Yoshino.peak.bloom.date.adj

stan_data <- list(N       = length(n_beta),
                  M       = M,
                  n_beta  = n_beta,
                  a       = a,
                  b       = b,
                  alpha_0 = alpha_0)

fit_bayesvern <- sampling(model, data=stan_data, init_r = 10, seed = 1)
fit_bayesvern
}

##Plots
M_bins <- cut(M,
              quantile(stan_data$M, seq(0, 1, length.out = 6)),
              include.lowest = TRUE)

sapply(levels(M_bins), function(bin) sd(n_beta[M_bins == bin]))

qplot(M_bins, n_beta, geom = "boxplot") +
  labs(x = "Average Temperature FEB - MAY",
       y = "Cherry Blossom Peak Bloom Date") +
  theme_bw()

# plot(Yoshino.peak.bloom.date~Year, data=cherry, type="l", col="hotpink", lwd=2)

## Eurodata
bp <- read.csv("..//pep_analyses/output/betpen_decsens_1950-2000.csv")
table(bp$siteslist)

bps1 <- subset(bp, siteslist==1)
bps2 <- subset(bp, siteslist==2)
bps3 <- subset(bp, siteslist==3) # picking a dramatic one, but other sites work too

## Plot of lm
bps3pre <- subset(bps3, year> 1949 & year<1971)
bps3post <- subset(bps3, year> 1989 & year<2011)

par(mfrow=c(1, 2))
plot(lo~mat, data=bps3pre)
abline(lm(lo~mat, data=bps3pre))
plot(lo~mat, data=bps3post)
abline(lm(lo~mat, data=bps3post))

## Stan

n_beta <- bps3$lo
M <- (bps3$mat * (9/5))+32

stan_data <- list(N       = length(n_beta),
                  M       = M,
                  n_beta  = n_beta,
                  a       = a,
                  b       = b,
                  alpha_0 = alpha_0)

##Plots
M_bins <- cut(M,
              quantile(stan_data$M, seq(0, 1, length.out = 6)),
              include.lowest = TRUE)

sapply(levels(M_bins), function(bin) sd(n_beta[M_bins == bin]))

qplot(M_bins, n_beta, geom = "boxplot") +
  labs(x = "Average Temperature March to May",
       y = "Betula pendula leafout date") +
  theme_bw()

fit_bayes <- sampling(model, data=stan_data, init_r = 10, seed = 1)
fit_bayes

