data {
  int<lower=0> N;                       // number of data points
  int<lower=0> Ng;                      // number of samples
  int<lower=0> Ns;                      // number of subjects
  int<lower=0> Nc;                      // number of conditions
  int<lower=0> Gcount [Ng];             // BCRs per sample
  real <lower=0> Y [N];                 // CDR3 lengths
  int <lower=0> G [N];                  // sample IDs
  int <lower=0> C [Ns];                 // condition IDs
  int <lower=0> S [Ng];                 // subject IDs
  real <lower=0> prior_mean;            // prior mean in condition (of normal pdf)
  real <lower=0> prior_sigma;           // prior sigma in condition (of normal pdf)
}


parameters {
  vector <lower=0> [Nc] mu_condition;
  vector <lower=0> [Ng] sigma_sample;
  vector <lower=0> [Ns] sigma_subject;
  vector <lower=0> [Nc] sigma_condition;
  vector [Ng] z_sample;
  vector [Ns] z_subject;
}


transformed parameters {
  vector <lower=0> [Ng] mu_sample;
  vector <lower=0> [Nc] mu_subject;

  for(g in 1:Ng) {
    mu_sample[g] = mu_subject[S[g]] + sigma_subject[S[g]]*z_sample[g];
  }

  for(s in 1:Ns) {
    mu_subject[s] = mu_condition[C[s]] + sigma_condition[C[s]]*z_subject[s];
  }
}

model {
  int pos;

  pos = 1;
  for(g in 1:Ng) {
    segment(Y, pos, Gcount[g]) ~ normal(mu_sample[g], sigma_sample[g]);
    pos = pos + Gcount[g];
  }

  mu_condition ~ normal(prior_mean, prior_sigma);

  sigma_sample ~ cauchy(0, 1);
  sigma_subject ~ cauchy(0, 1);
  sigma_condition ~ cauchy(0, 1);

  z_sample ~ std_normal();
  z_subject ~ std_normal();
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

