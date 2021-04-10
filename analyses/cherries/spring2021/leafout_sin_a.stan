/* Notes by Lizzie (9 April 2021) */
/* But code by J Auerbach */

/* Similar to leafout_sin_p.stan but ...
Here we assume plants start accumulating on a (say January 1) and don't try to caculate it (known a)
In the leafout_sin_p.stan we add p and try to estimate it */

functions {
  real temp_sum(int d, real alpha_2) {
    return
       sin(2 * pi() * alpha_2 / 365) *
         (sin(pi() / 365) - sin(pi() * (2 * d + 1) / 365))  + 
       cos(2 * pi() * alpha_2 / 365) * 
         (cos(pi() / 365) - cos(pi() * (2 * d + 1) / 365));
  }
  real n_given_m(int n, real m, real alpha_0, real alpha_1, real alpha_2, 
                 real sigma, int a, int b, real s) {
    return
       Phi(
           (b * m - s - (
             alpha_0 * (b - n + a) + alpha_1 / (2 * sin(pi() / 365)) * 
             (temp_sum(b, alpha_2) - temp_sum(n, alpha_2) + temp_sum(a, alpha_2))
            ) 
           ) / (sigma * sqrt(b + a - n))
          );
    }
  real log_lik(int n, real m, real alpha_0, real alpha_1, real alpha_2,
                real sigma, int a, int b, real s) {
    return 
      log(n_given_m(n, m, alpha_0, alpha_1, alpha_2, sigma, a, b, s) - 
          n_given_m(n - 1, m, alpha_0, alpha_1, alpha_2, sigma, a, b, s)) +
      normal_lpdf(m | alpha_0 + 
                      alpha_1 * temp_sum(b, alpha_2) / (2 * b * sin(pi() / 365)), 
                      sigma / sqrt(b));
    }
}

data {
  int<lower=0> N;           // number of experiments i = 1, ..., N
  int a;                    // start of plant accumualation 
  int b;                    // end of average temperature of window
  vector[N] M;              // mean temperature from 0 to b for each experiment
  int n_s[N];               // leaf out day for each experiment
}

parameters {
  real<lower = 0> alpha_0;
  real<lower = 0, upper = alpha_0/2> alpha_1;
  real<lower = 90, upper = 120>      alpha_2;
  real<lower = 0> s;
  real<lower = 1> sigma;
}

model {
  for(n in 1:N) {
    target += log_lik(n_s[n], M[n], alpha_0, alpha_1, alpha_2, sigma, a, b, s);
  }
}
