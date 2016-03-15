// LUX
data {
  real g_a;
  real g_b;

  int<lower=1> R;
  int<lower=1> N;

  real<lower=0,upper=1> bsEff[R];
  real<lower=0,upper=1> bsBEff[R];
  real<lower=0,upper=1> oxEff[R];
  real<lower=0,upper=1> seqErr[R];

  int<lower=0> bsC[N,R];
  int<lower=0> bsTot[N,R];
  int<lower=0> oxC[N,R];
  int<lower=0> oxTot[N,R];

  vector<lower=0>[3] alpha[N];
}
parameters {
  real<lower=0> g[N];

  simplex[3] mu[N];
  simplex[3] theta[N,R];
}
transformed parameters {
  vector<lower=0>[3] mug[N];

  for (n in 1:N) {
    mug[n] <- g[n]*mu[n]+1;
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
  }

  for (n in 1:N) {
    mu[n] ~ dirichlet(alpha[n]);
  }

  g ~ gamma(g_a,g_b);
}
