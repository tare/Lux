// LUX
data {
  real g_a;
  real g_b;

  real bsEff;
  real bsBEff;
  real oxEff;
  real labEff;
  real proEff;
  real seqErr;

  int<lower=1> R;
  int<lower=1> N;

  int<lower=0> bsC[R,N];
  int<lower=0> bsTot[R,N];
  int<lower=0> tabC[R,N];
  int<lower=0> tabTot[R,N];
  int<lower=0> fcabC[R,N];
  int<lower=0> fcabTot[R,N];

  vector<lower=0>[4] alpha[N];
}
parameters {
  real<lower=0> g[N];

  simplex[4] mu[N];

  simplex[4] theta[N,R];
}
transformed parameters {
  vector<lower=0>[4] mug[N];

  for (n in 1:N) {
    mug[n] <- g[n]*mu[n]+1;
  }
}
model {
  real p_BS[R,N];
  real p_TAB[R,N];
  real p_fCAB[R,N];

  for (r in 1:R) {
    for (n in 1:N) {
      theta[n,r] ~ dirichlet(mug[n]);
      p_BS[r,n] <- theta[n,r,1]*((1.0 - bsEff)*(1.0 - seqErr) + seqErr * bsEff) +
                   theta[n,r,2]*((1.0 - bsBEff)*(1.0 - seqErr) + seqErr * bsBEff) +
                   theta[n,r,3]*((1.0 - bsBEff)*(1.0 - seqErr) + seqErr * bsBEff) +
                   theta[n,r,4]*((1.0 - bsEff)*(1.0 - seqErr) + seqErr * bsEff);

      p_TAB[r,n] <- theta[n,r,1]*((1.0 - bsEff)*(1.0 - seqErr) + seqErr * bsEff) +
                    theta[n,r,2]*(oxEff * ((1.0-bsEff) * (1.0-seqErr) + bsEff * seqErr) + (1.0 - oxEff) * ((1.0 - bsBEff) * (1.0 - seqErr) + bsBEff * seqErr)) +
                    theta[n,r,3]*(labEff * ((1.0 - bsBEff) * (1.0 - seqErr) + bsBEff * seqErr) + (1.0 - labEff) * (oxEff * ((1.0 - bsEff) * (1.0 - seqErr) + bsEff * seqErr) + (1.0 - oxEff) * ((1.0 - bsBEff) * (1.0 - seqErr) + bsBEff * seqErr))) +
                    theta[n,r,4]*(oxEff * ((1.0 - bsEff) * (1.0 - seqErr) + bsEff * seqErr) + (1.0 - oxEff) * ((1.0 - bsEff) * (1.0 - seqErr) + bsEff * seqErr));

      p_fCAB[r,n] <- theta[n,r,1]*((1.0 - bsEff)*(1.0 - seqErr) + seqErr * bsEff) +
                     theta[n,r,2]*((1.0 - bsBEff)*(1.0 - seqErr) + seqErr * bsBEff) +
                     theta[n,r,3]*((1.0 - bsBEff)*(1.0 - seqErr) + seqErr * bsBEff) +
                     theta[n,r,4]*(proEff * ((1.0 - bsBEff) * (1.0 - seqErr) + bsBEff * seqErr) + (1.0 - proEff) * ((1.0 - bsEff) * (1.0 - seqErr) + bsEff * seqErr));
    }
  }

  for (r in 1:R) {
    bsC[r] ~ binomial(bsTot[r], p_BS[r]);
    tabC[r] ~ binomial(tabTot[r], p_TAB[r]);
    fcabC[r] ~ binomial(fcabTot[r], p_fCAB[r]);
  }

  for (n in 1:N) {
    mu[n] ~ dirichlet(alpha[n]);
  }

  g ~ gamma(g_a,g_b);
}
