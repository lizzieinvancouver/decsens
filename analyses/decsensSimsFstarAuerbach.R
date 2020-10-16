## Wrote 15 Oct 2020 ##
## From J Auerbach's 12 Oct 2020 email ##


#Plot 1 (model from section 1.1)

# Example plot for Section 2.1
thresholds <- seq(100, 400, length.out = 10)
cov_experiment <- var_experiment <- numeric(length(thresholds))
for(experiment in seq_along(cov_experiment)) {
  mean_temp <- leaf_out <- numeric(1e5)
  for(year in seq_along(mean_temp)) {
    temp <- rnorm(100, 5)
    leaf_out[year] <- which.min(!cumsum(temp) > thresholds[experiment])
    mean_temp[year] <- mean(temp[1:leaf_out[year]])
  }
  cov_experiment[experiment] <- cov(mean_temp, leaf_out)
  var_experiment[experiment] <- var(mean_temp)
}

par(mar = c(5, 5, 3, 5))
plot(thresholds, cov_experiment, axes=F,
     ylim=range(cov_experiment),
     xlab="threshold", ylab = "covariance", type="l",
     col="blue", main="",xlim=range(thresholds))
axis(side = 1)
axis(side = 2)
par(new=T)
plot(thresholds, var_experiment, axes=F,
     ylim=range(var_experiment), type="l", lty = 2,
     xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     col="red", main="",xlim=range(thresholds))
axis(side = 4)
mtext("variance", side = 4, line = 3)
legend("topleft", c("covariance", "variance"),
       col = c("blue", "red"), lty = c(1, 2))

#Plot 2 (model from section 1.2)


thresholds <- seq(100, 400, length.out = 10)
cov_experiment <- var_experiment <- numeric(length(thresholds))
for(experiment in seq_along(cov_experiment)) {
  mean_temp <- leaf_out <- numeric(1e5)
  for(year in seq_along(mean_temp)) {
    temp <- rnorm(100, 5)
    leaf_out[year] <- which.min(!cumsum(temp) > thresholds[experiment])
    mean_temp[year] <- mean(temp)
  }
  cov_experiment[experiment] <- cov(mean_temp, leaf_out)
  var_experiment[experiment] <- var(mean_temp)
}

par(mar = c(5, 5, 3, 5))
plot(thresholds, cov_experiment, axes=F,
     ylim=range(cov_experiment),
     xlab="threshold", ylab = "covariance", type="l",
     col="blue", main="",xlim=range(thresholds))
axis(side = 1)
axis(side = 2)
par(new=T)
plot(thresholds, var_experiment, axes=F,
     ylim=range(var_experiment), type="l", lty = 2,
     xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     col="red", main="",xlim=range(thresholds))
axis(side = 4)
mtext("variance", side = 4, line = 3)
legend("topleft", c("covariance", "variance"),
       col = c("blue", "red"), lty = c(1, 2))
