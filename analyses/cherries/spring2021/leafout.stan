functions {
  real n_given_m(int n, real m, real alpha_0, real alpha_1, 
                 real sigma, int b, real p, real s) {
    real n_and_a_given_m = 0;
    for(a in 1:b) {
      n_and_a_given_m += 
        Phi(
            (b * m - s - (alpha_0 * (b - n + a) + 
               alpha_1 / 2 * ( b * (b + 1) - n * (n + 1) + a * (a + 1))) 
            ) / (sigma * sqrt(b + a - n))
           ) * (p * (1 - p) ^ (a - 1)) / (1 - (1 - p) ^ b);
    }
    return n_and_a_given_m;
  }
  real log_lik(int n, real m, real alpha_0, real alpha_1,
                real sigma, int b, real p, real s) {
    return 
      log(n_given_m(n, m, alpha_0, alpha_1, sigma, b, p, s) - 
          n_given_m(n - 1, m, alpha_0, alpha_1, sigma, b, p, s)) +
      normal_lpdf(m | alpha_0 + alpha_1 / 2 * (b + 1), 
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
  real<lower = 0> alpha_1;
  real<lower = 0> s;
  real<lower = 1> sigma;
  real<lower = 0, upper = 1> p;
}

model {
  for(n in 1:N) {
    target += log_lik(n_s[n], M[n], alpha_0, alpha_1, sigma, b, p, s);
  }
}
