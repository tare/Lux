data {
  real mu_mu_bsEff;
  real sigma_mu_bsEff;
  
  real mu_sigma_bsEff;
  real sigma_sigma_bsEff;

  real mu_mu_bsBEff;
  real sigma_mu_bsBEff;

  real mu_sigma_bsBEff;
  real sigma_sigma_bsBEff;

  real mu_mu_oxEff;
  real sigma_mu_oxEff;

  real mu_sigma_oxEff;
  real sigma_sigma_oxEff;

  real mu_mu_seqErr;
  real sigma_mu_seqErr;

  real mu_sigma_seqErr;
  real sigma_sigma_seqErr;

  real g_a;
  real g_b;

  int<lower=1> R;
  int<lower=1> N;
  int<lower=1> N_control;

  int<lower=0> bsC[N,R];
  int<lower=0> bsTot[N,R];
  int<lower=0> oxC[N,R];
  int<lower=0> oxTot[N,R];

  int<lower=0> bsC_control[N_control,R];
  int<lower=0> bsTot_control[N_control,R];
  int<lower=0> oxC_control[N_control,R];
  int<lower=0> oxTot_control[N_control,R];

  vector<lower=0>[3] alpha[N];

  vector<lower=0>[3] alpha_control[N_control];
}
parameters {
  real mu_bsEff;
  real<lower=0> sigma_bsEff;
  real raw_bsEff[R];

  real mu_bsBEff;
  real<lower=0> sigma_bsBEff;
  real raw_bsBEff[R];

  real mu_oxEff;
  real<lower=0> sigma_oxEff;
  real raw_oxEff[R];

  real mu_seqErr;
  real<lower=0> sigma_seqErr;
  real raw_seqErr[R];

  real<lower=0> g[N];

  simplex[3] mu[N];
  simplex[3] theta[N,R];

  simplex[3] theta_control[N_control,R];
}
transformed parameters {
  real<lower=0,upper=1> bsEff[R];
  real<lower=0,upper=1> bsBEff[R];
  real<lower=0,upper=1> oxEff[R];
  real<lower=0,upper=1> seqErr[R];

  vector<lower=0>[3] mug[N];

  for (n in 1:N) {
    mug[n] <- g[n]*mu[n]+1;
  }

  for (r in 1:R) {
    bsEff[r]  <- inv_logit(mu_bsEff  + sigma_bsEff  * raw_bsEff[r]);
    bsBEff[r] <- inv_logit(mu_bsBEff + sigma_bsBEff * raw_bsBEff[r]);
    oxEff[r]  <- inv_logit(mu_oxEff  + sigma_oxEff  * raw_oxEff[r]);
    seqErr[r] <- inv_logit(mu_seqErr + sigma_seqErr * raw_seqErr[r]);
  }
}
model {
  for (r in 1:R) {
    for (n in 1:N) {
      theta[n,r] ~ dirichlet(mug[n]);
      bsC[n,r] ~ binomial(bsTot[n,r],
              theta[n,r,1]*((1.0 - seqErr[r])*(1.0 - bsEff[r]) + seqErr[r] * bsEff[r]) +
              theta[n,r,2]*((1.0 - bsBEff[r])*(1.0 - seqErr[r]) + seqErr[r] * bsBEff[r]) +
              theta[n,r,3]*((1.0 - bsBEff[r])*(1.0 - seqErr[r]) + seqErr[r] * bsBEff[r]));
      oxC[n,r] ~ binomial(oxTot[n,r],
              theta[n,r,1]*((1.0 - seqErr[r])*(1.0 - bsEff[r]) + seqErr[r] * bsEff[r]) +
              theta[n,r,2]*((1.0 - bsBEff[r])*(1.0 - seqErr[r]) + seqErr[r] * bsBEff[r]) +
              theta[n,r,3]*(oxEff[r]*((1.0 - seqErr[r])*(1.0 - bsEff[r]) + seqErr[r] * bsEff[r])+(1.0 - oxEff[r])*((1.0 - bsBEff[r])*(1.0 - seqErr[r]) + bsBEff[r] * seqErr[r])));
    }
    for (n in 1:N_control) {
      theta_control[n,r] ~ dirichlet(alpha_control[n]);
      bsC_control[n,r] ~ binomial(bsTot_control[n,r],
              theta_control[n,r,1]*((1.0 - seqErr[r])*(1.0 - bsEff[r]) + seqErr[r] * bsEff[r]) +
              theta_control[n,r,2]*((1.0 - bsBEff[r])*(1.0 - seqErr[r]) + seqErr[r] * bsBEff[r]) +
              theta_control[n,r,3]*((1.0 - bsBEff[r])*(1.0 - seqErr[r]) + seqErr[r] * bsBEff[r]));
      oxC_control[n,r] ~ binomial(oxTot_control[n,r],
              theta_control[n,r,1]*((1.0 - seqErr[r])*(1.0 - bsEff[r]) + seqErr[r] * bsEff[r]) +
              theta_control[n,r,2]*((1.0 - bsBEff[r])*(1.0 - seqErr[r]) + seqErr[r] * bsBEff[r]) +
              theta_control[n,r,3]*(oxEff[r]*((1.0 - seqErr[r])*(1.0 - bsEff[r]) + seqErr[r] * bsEff[r])+(1.0 - oxEff[r])*((1.0 - bsBEff[r])*(1.0 - seqErr[r]) + bsBEff[r] * seqErr[r])));
    }
  }

  for (n in 1:N) {
    mu[n] ~ dirichlet(alpha[n]);
  }

  g ~ gamma(g_a,g_b);

  mu_bsEff ~ normal(mu_mu_bsEff, sigma_mu_bsEff);
  sigma_bsEff ~ lognormal(mu_sigma_bsEff, sigma_sigma_bsEff);
  raw_bsEff ~ normal(0,1);

  mu_bsBEff ~ normal(mu_mu_bsBEff, sigma_mu_bsBEff);
  sigma_bsBEff ~ lognormal(mu_sigma_bsBEff, sigma_sigma_bsBEff);
  raw_bsBEff ~ normal(0,1);

  mu_oxEff ~ normal(mu_mu_oxEff, sigma_mu_oxEff);
  sigma_oxEff ~ lognormal(mu_sigma_oxEff, sigma_sigma_oxEff);
  raw_oxEff ~ normal(0,1);

  mu_seqErr ~ normal(mu_mu_seqErr, sigma_mu_seqErr);
  sigma_seqErr ~ lognormal(mu_sigma_seqErr, sigma_sigma_seqErr);
  raw_seqErr ~ normal(0,1);
}
