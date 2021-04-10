functions {
  real n_beta_given_m_cdf(int n, real m, real alpha_0, real alpha_1, 
                          real sigma, int b, int a, real beta) {
  return Phi(
   (
    ((b - a) * m - beta) - 
     (alpha_0 * (b - a - n) + alpha_1/2 * ((b - n) * (b + n + 1) - a * (a + 1)))
   ) / (sigma * sqrt(b + a - n))
   );
  }
}

data {
  int<lower=0> N;           // number of experiments n = 1, ..., N
  int<lower=0> S;           // number of sites
  int<lower=0> P;           // number of covariates
  int site[N];              // site of experiment n
  vector[N] M;              // mean temperature in experiment n
  int n_beta[N];            // leaf out day in experiment n
  matrix[N, P] X;           // matrix of P scaled covariates for experiment n
  int a;                    // beginning of average temperature window 
  int b;                    // end of average temperature of window
}

parameters {
  real<lower = 0> sigma[S];           // standard deviation of daily temperature
  vector[P-1] gamma_alpha_0;          // regression coefficients: baseline temp
  vector[P]   gamma_alpha_1;          //                          daily temp inc
  vector[P]   gamma_beta;             //                          temp threshold
  real<lower = 0> tau_gamma_alpha_0;  // standard deviation of coef: baseline
  real<lower = 0> tau_gamma_alpha_1;  //                             daily incr 
  real<lower = 0> tau_gamma_beta;     //                             temp thresh
}

transformed parameters {
  vector[N] alpha_0 = 32.0 * exp(block(X, 1, 2, N, P-1) * gamma_alpha_0);
  vector[N] alpha_1 = exp(X * gamma_alpha_1);
  vector[N] beta    = exp(X * gamma_beta);
}

model {
  for(n in 1:N) {
    target += 
//log likelihood of n_beta | M    
    log(n_beta_given_m_cdf(n_beta[n], M[n], 
                           alpha_0[n], alpha_1[n], 
                           sigma[site[n]], b, a, beta[n]) -
                  n_beta_given_m_cdf(n_beta[n] - 1, M[n], 
                           alpha_0[n], alpha_1[n], 
                           sigma[site[n]], b, a, beta[n])) +      
//log likelihood of M
               normal_lpdf(M[n] | 
                           alpha_0[n] + alpha_1[n]/2 * (b + a + 1), 
                           sigma[site[n]] / sqrt(b-a));
  }
//regularization of coefficients (excluding intercept)
  for(p in 2:P) {
    gamma_alpha_0[p-1] ~ normal(0, tau_gamma_alpha_0);
    gamma_alpha_1[p]   ~ normal(0, tau_gamma_alpha_1);
    gamma_beta[p]      ~ normal(0, tau_gamma_beta);
  }
}
