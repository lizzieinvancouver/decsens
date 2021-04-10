/* Notes by Lizzie (9 April 2021) */
/* But code by J Auerbach */


/* Similar to leafout_sin_a.stan but ...
Here we assume we don't know when plants start accumulating, so we p and try to estimate it */

functions {
  real temp_sum(int d, real alpha_2) {
    return
       sin(2 * pi() * alpha_2 / 365) *
         (sin(pi() / 365) - sin(pi() * (2 * d + 1) / 365))  + 
       cos(2 * pi() * alpha_2 / 365) * 
         (cos(pi() / 365) - cos(pi() * (2 * d + 1) / 365));
  }
  real n_given_m(int n, real m, real alpha_0, real alpha_1, real alpha_2, 
                 real sigma, int b, real p, real s) {
    real n_and_a_given_m = 0;
/* below chunk is marginalizing out since a is discrete parameter (days) 
 this is also a spot where the code does a LOT of summations */ 
    for(a in 1:b) { 
      n_and_a_given_m += 
       Phi(
           (b * m - s - (
             alpha_0 * (b - n + a) + alpha_1 / (2 * sin(pi() / 365)) * 
             (temp_sum(b, alpha_2) - temp_sum(n, alpha_2) + temp_sum(a, alpha_2))
            ) 
           ) / (sigma * sqrt(b + a - n))
          ) * (p * (1 - p) ^ (a - 1)) / (1 - (1 - p) ^ b);
    }
    return n_and_a_given_m;
  }
  real log_lik(int n, real m, real alpha_0, real alpha_1, real alpha_2,
                real sigma, int b, real p, real s) {
    return 
      log(n_given_m(n, m, alpha_0, alpha_1, alpha_2, sigma, b, p, s) - 
          n_given_m(n - 1, m, alpha_0, alpha_1, alpha_2, sigma, b, p, s)) +
      normal_lpdf(m | alpha_0 + 
                      alpha_1 * temp_sum(b, alpha_2) / (2 * b * sin(pi() / 365)), 
                      sigma / sqrt(b));
    }
}

data {
  int<lower=0> N;           // number of experiments i = 1, ..., N
  int b;                    // end of average temperature of window
  vector[N] M;              // mean temperature from 0 to b for each experiment
  int n_s[N];               // leaf out day for each experiment
}

parameters {
  real<lower = 0> alpha_0;
  real<lower = 0, upper = alpha_0/2> alpha_1; // amplitude of temperature curve
  real<lower = 90, upper = 120>      alpha_2; // phase (moves it over so temperature go up at a reasonable day of year
  real<lower = 0> s;
  real<lower = 1> sigma; 
  real<lower = 0, upper = 1> p; // a coin toss ... when plants start to accumulate temperatures
}

model {
  for(n in 1:N) {
    target += log_lik(n_s[n], M[n], alpha_0, alpha_1, alpha_2, sigma, b, p, s);
  }
}
