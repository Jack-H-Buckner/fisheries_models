// Bayesian schnutes method
data {
  // data
  int<lower=0> N;
  vector[N] y;
  vector[N] x1;
  vector[N] x2;
  
  // priors
  real alpha_r;
  real beta_r;
  real alpha_q;
  real beta_q;
  real alpha_K;
  real beta_K;
}

// regression parameters
// y = a + bx1 + cx2
parameters {
  real<lower=0> a;
  real<upper=0> b;
  real<upper=0> c;
  real<lower=0> sigma;
}
// transform regression parameters to 
// biological parameters
transformed parameters {
  real<lower=0> r;
  real<lower=0> K;
  real<lower=0> q;
  r = a;
  q = -1*c;
  K = a/(c*b);
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // priors
  target += gamma_lpdf(r|alpha_r,beta_r);
  target += gamma_lpdf(q|alpha_q,beta_q);
  target += gamma_lpdf(K|alpha_K,beta_K);
  // likelihood
  y ~ normal(a + b*x1 + c*x2, sigma);
}





