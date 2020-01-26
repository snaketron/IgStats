data {
  int<lower=0> N;                       // number of data points
  int<lower=0> Ng;                      // number of all samples
  int<lower=0> Nc;                      // number of conditions
  int<lower=0> Gcount [Ng];             // number of groups
  real Y [N];                           // CDR3 charge
  int <lower=0> C [Ng];                 // condition ID
  int <lower=0> G [N];                  // sample IDs
  real prior_mean;                      // prior mean in condition (of normal pdf)
  real <lower=0> prior_sigma;           // prior sigma in condition (of normal pdf)
}

parameters {
  vector <lower=0> [Ng] sigma_sample;
  vector [Nc] mu_condition;
  vector <lower=0> [Nc] sigma_condition;
  vector [Ng] z_sample;
}

transformed parameters {
  vector [Ng] mu_sample;
  for(g in 1:Ng) {
    mu_sample[g] = mu_condition[C[g]] + sigma_condition[C[g]]*z_sample[g];
  }
}

model {
  int pos;

  pos = 1;
  for(g in 1:Ng) {
    segment(Y, pos, Gcount[g]) ~ normal(mu_sample[g], sigma_sample[g]);
    pos = pos + Gcount[g];
  }

  for(c in 1:Nc) {
    mu_condition[c] ~ normal(prior_mean, prior_sigma);
  }

  sigma_condition ~ cauchy(0, 1);
  sigma_sample ~ cauchy(0, 1);
  z_sample ~ std_normal();
}


generated quantities {
  vector [Ng] Yhat_sample;
  real log_lik [N];

  for(g in 1:Ng) {
    Yhat_sample[g] = normal_rng(mu_sample[g], sigma_sample[g]);
  }

  for(i in 1:N) {
    log_lik[i] = normal_lpdf(Y[i] | mu_sample[G[i]], sigma_sample[G[i]]);
  }
}

