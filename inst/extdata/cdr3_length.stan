data {
  int<lower=0> N;                       // number of data points
  int<lower=0> Ng;                      // number of all samples
  int<lower=0> Nc;                      // number of conditions
  int<lower=0> Gcount [Ng];             // number of groups
  real <lower=0> Y [N];                 // CDR3 length
  int <lower=0> C [Ng];                 // condition ID
  real <lower=0> prior_mean_condition;  // prior mean in condition (of normal pdf)
  real <lower=0> prior_sigma_condition; // prior sigma in condition (of normal pdf)
}


parameters {
  vector <lower=0> [Ng] mu_sample;
  vector <lower=0> [Ng] sigma_sample;
  vector <lower=0> [Nc] mu_condition;
  vector <lower=0> [Nc] sigma_condition;
}


model {
  int pos;

  pos = 1;
  for(g in 1:Ng) {
    segment(Y, pos, Gcount[g]) ~ normal(mu_sample[g], sigma_sample[g]);
    pos = pos + Gcount[g];
    mu_sample[g] ~ normal(mu_condition[C[g]], sigma_condition[C[g]]);
  }

  for(c in 1:Nc) {
    mu_condition[c] ~ normal(prior_mean_condition, prior_sigma_condition);
  }

  sigma_condition ~ cauchy(0, 1);
  sigma_sample ~ cauchy(0, 1);
}


generated quantities {
  vector [Ng] Yhat_sample;
  real rng;

  for(g in 1:Ng) {
    rng = -1;
    while(rng < 0) {
      rng = normal_rng(mu_sample[g], sigma_sample[g]);
    }
    Yhat_sample[g] = rng;
  }
}
