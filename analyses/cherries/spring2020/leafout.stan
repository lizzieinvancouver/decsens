functions {
  real n_beta_given_m_cdf(int n, real m, 
                   real alpha_0, real alpha_1, real sigma, 
                   int b, int a, real beta) {
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
  vector[N] M;              // mean temperature in year n
  int n_beta[N];            // leaf out day in year n
  int a;                    // beginning of average temperature window 
  int b;                    // end of average temperature of window
  real alpha_0;             // baseline temperature on day 0
}

parameters {
  real<lower = 0> alpha_1;
  real<lower = 0, upper = 1> beta_unit_scale;
  real<lower = 1> sigma;
}

transformed parameters {
  real beta = (b * max(M)) * beta_unit_scale + 
               (alpha_0 * min(n_beta)) * (1 - beta_unit_scale);
//n.b. this is equivalentish to 
 // real<lower = b * max(M), upper= alpha_0 * min(n_beta)> beta
}

model {
  for(n in 1:N) {
    target += 
//log likelihood of n_beta | M    
    log(n_beta_given_m_cdf(n_beta[n], M[n], 
                           alpha_0, alpha_1, sigma, 
                           b, a, beta) -
                  n_beta_given_m_cdf(n_beta[n] - 1, M[n], 
                           alpha_0, alpha_1, sigma, 
                           b, a, beta)) +      
//log likelihood of M
               normal_lpdf(M[n] | 
                           alpha_0 + alpha_1/2 * (b + a + 1), 
                           sigma / sqrt(b-a));
  }
}
