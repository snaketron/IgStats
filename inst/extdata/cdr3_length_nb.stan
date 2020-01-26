data {
  int<lower=0> N;                       // number of data points
  int<lower=0> Ng;                      // number of all samples
  int<lower=0> Nc;                      // number of conditions
  int<lower=0> Gcount [Ng];             // number of groups
  int <lower=0> Y [N];                  // CDR3 length
  int <lower=0> C [Ng];                 // condition ID
  int <lower=0> G [N];                  // sample IDs
  real <lower=0> prior_mean_condition;  // prior mean in condition (of normal pdf)
  real <lower=0> prior_sigma_condition; // prior sigma in condition (of normal pdf)
}


parameters {
  vector <lower=0> [Ng] phi_sample;
  vector <lower=0> [Ng] tau_sample;
  vector [Ng] z_sample;
  vector <lower=0> [Nc] mu_condition;
  vector <lower=0> [Nc] sigma_condition;
}


transformed parameters {
  vector <lower=0> [Ng] mu_sample;

  for(g in 1:Ng) {
    mu_sample[g] = mu_condition[C[g]] + sigma_condition[C[g]]*z_sample[g];
  }
}


model {
  int pos;

  pos = 1;
  for(g in 1:Ng) {
    segment(Y, pos, Gcount[g]) ~ neg_binomial_2(mu_sample[g], phi_sample[g]);
    pos = pos + Gcount[g];
  }

  phi_sample ~ exponential(tau_sample);
  tau_sample ~ gamma(2.0, 0.1);
  mu_condition ~ normal(prior_mean_condition, prior_sigma_condition);
  sigma_condition ~ cauchy(0, 1);
  z_sample ~ std_normal();
}


generated quantities {
  vector [Ng] Yhat_sample;
  real log_lik [N];

  for(g in 1:Ng) {
    Yhat_sample[g] = neg_binomial_2_rng(mu_sample[g], phi_sample[g]);
  }

  for(i in 1:N) {
    log_lik[i] = neg_binomial_2_lpmf(Y[i] | mu_sample[G[i]], phi_sample[G[i]]);
  }
}
