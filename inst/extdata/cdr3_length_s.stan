data {
  int<lower=0> N;                       // number of data points
  int<lower=0> Ng;                      // number of samples
  int<lower=0> Gcount [Ng];             // number of groups
  real <lower=0> Y [N];                 // CDR3 length
  int <lower=0> G [N];                  // sample IDs
  real <lower=0> prior_mean;            // prior mean in samples (of normal pdf)
  real <lower=0> prior_sigma;           // prior sigma in samples (of normal pdf)
}


parameters {
  vector <lower=0> [Ng] mu_sample;
  vector <lower=0> [Ng] sigma_sample;
}


model {
  int pos;

  pos = 1;
  for(g in 1:Ng) {
    segment(Y, pos, Gcount[g]) ~ normal(mu_sample[g], sigma_sample[g]);
    pos = pos + Gcount[g];
  }

  mu_sample ~ normal(prior_mean, prior_sigma);
  sigma_sample ~ cauchy(0, 1);
}


generated quantities {
  vector [Ng] Yhat_sample;
  // real log_lik [N];
  real rng;

  for(g in 1:Ng) {
    rng = -1;
    while(rng < 0) {
      rng = normal_rng(mu_sample[g], sigma_sample[g]);
    }
    Yhat_sample[g] = rng;
  }

  // for(i in 1:N) {
  //   log_lik[i] = normal_lpdf(Y[i] | mu_sample[G[i]], sigma_sample[G[i]]);
  // }
}

