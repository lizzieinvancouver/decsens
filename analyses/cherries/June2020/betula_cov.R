library("tidyverse")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Set temperature window
a <- 32            # Feb 1st, beginning of window for average temperature
b <- 121           # May 31st, end of window for average temperature

#################
# Read Eurodata #
#################
bp <- read.csv("betpen_decsens_1950-2000.csv")

#############
# Model Fit #
#############
bp_stan <- subset(bp, siteslist <= 17 & year > 1951)
n_beta <- bp_stan$lo
M      <- (bp_stan$mat * (9/5))+32
year   <- as.numeric(factor(bp_stan$year))
site   <- as.numeric(factor(bp_stan$siteslist))
X      <- scale(as.matrix(model.matrix(lo ~ factor(year) + 
                                            factor(bp_stan$siteslist) +
                                            gdd + chillutah, 
                                       bp_stan)))
X[ , 1] <- 1

stan_data <- list(N       = length(n_beta),
                  S       = length(unique(site)),
                  P       = ncol(X),
                  site    = site,
                  M       = M,
                  n_beta  = n_beta,
                  X       =  X,
                  a       = a,
                  b       = b)

init_cov <- function () 
  list(alpha_0         = rep(32, stan_data$N),
       alpha_1         = rep(.5, stan_data$N),
       beta            = rep(b * max(M), stan_data$N),
       sigma           = rep(100, stan_data$S),
       gamma_alpha_0     = rep(0, stan_data$P-1),
       gamma_alpha_1     = c(log(.5), rep(0, stan_data$P-1)),
       gamma_beta        = c(log(b * max(M)), rep(0, stan_data$P-1)),
       tau_gamma_alpha_0 = 100,
       tau_gamma_alpha_1 = 100,
       tau_gamma_beta    = 100
  )

model_cov <- stan_model("betula_cov.stan")
fit_bayes_cov <- sampling(model_cov, data = stan_data, init = init_cov, 
                          control = list(adapt_delta = 0.9))

##############
# Model Plot #
##############
#evidence of declining sensitivity after adjusting for warming temperatures
tibble(
    `alpha[0]` = colMeans(rstan::extract(fit_bayes_cov, "alpha_0")[[1]]),
    `alpha[1]` = colMeans(rstan::extract(fit_bayes_cov, "alpha_1")[[1]]),
    beta       = colMeans(rstan::extract(fit_bayes_cov, "beta")[[1]]),
    year       = year) %>%
  gather(key = "key", value = "value", -year) %>%
  ggplot() +
    aes(x = year, y = value) +
    theme_bw() +
    geom_point() +
    geom_smooth(aes(color = key), method = "lm",
                formula = y ~ splines::ns(x, knots = c(20, 40))) +
    facet_grid(key ~ ., scales = "free", labeller = label_parsed, switch = "y") +
    theme(strip.text.y = element_text(angle=180),
          legend.position = "none") +
    scale_y_continuous(name = "", position = "right") +
    scale_x_continuous(labels = function(x) x + min(bp$year)) +
    labs(x = "", 
         title = "Warming (top/middle) and Sensitivity (bottom) for Betula pendula")

#regression coefficients
get_quantile <- function(var, param, probs) 
  if(param == "gamma_alpha_0") {
  exp(quantile(
    extract(fit_bayes_cov, param)[[1]][,which(colnames(X) == var)-1], 
    probs)) } else {
  exp(quantile(
    extract(fit_bayes_cov, param)[[1]][,which(colnames(X) == var)], 
    probs))
      }

tibble(
  probs = c(.025, .25, .5, .75, .975),
  `gdd_alpha[0]` = get_quantile("gdd", "gamma_alpha_0", probs),
  `gdd_alpha[1]` = get_quantile("gdd", "gamma_alpha_1", probs),
   gdd_beta      = get_quantile("gdd", "gamma_beta", probs),
  `chillutah_alpha[0]` = get_quantile("chillutah", "gamma_alpha_0", probs),
  `chillutah_alpha[1]` = get_quantile("chillutah", "gamma_alpha_1", probs),
   chillutah_beta      = get_quantile("chillutah", "gamma_beta", probs)
) %>%
  gather(key = key, value = value, -probs) %>%
  mutate(var = str_split(key, "_", simplify = TRUE)[, 1],
         param = str_split(key, "_", simplify = TRUE)[, 2]) %>%
  spread(key = probs, value = value) %>%
  ggplot() +
   aes(x = var, y = `0.5`) +
   theme_bw() +
   geom_linerange(aes(ymin = `0.25`, ymax = `0.75`), size = 2) +
   geom_linerange(aes(ymin = `0.025`, ymax = `0.975`))  +
   geom_point(color = "red") +
   facet_grid(param ~ ., scales = "free", labeller = label_parsed, switch = "y") +
   labs(x = "", y = "multiplicative effect: median, inner 50%, inner 95%", 
        title = "Warming (top/middle) and Sensitivity (bottom) for Betula pendula") +
   geom_hline(yintercept = 1, linetype = 2) +
   scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(strip.text.y = element_text(angle=180),
        legend.position = "none")