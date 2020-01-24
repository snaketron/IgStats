data {
  int<lower=0> N;                       // number of data points
  int<lower=0> Ng;                      // number of all samples
  int<lower=0> Nc;                      // number of conditions
  int<lower=0> Np;                      // number of productive statuses
  int<lower=0> Gcount [Ng];             // number of groups
  real Y [N];                           // CDR3 length
  int <lower=0> P [Ng];                 // productive or non-productive
  int <lower=0> C [Ng];                 // condition ID
  real prior_mean_condition;            // prior mean in condition (of normal pdf)
  real <lower=0> prior_sigma_condition; // prior sigma in condition (of normal pdf)
}


parameters {
  vector [Ng] mu_sample;
  vector <lower=0> [Ng] sigma_sample;
  vector [Nc] mu_condition [Np];
  vector <lower=0> [Nc] sigma_condition [Np];
}


model {
  int pos;
  
  pos = 1;
  for(g in 1:Ng) {
    segment(Y, pos, Gcount[g]) ~ normal(mu_sample[g], sigma_sample[g]);
    pos = pos + Gcount[g];
    mu_sample[g] ~ normal(mu_condition[P[g]][C[g]], sigma_condition[P[g]][C[g]]);
  }
  
  
  for(p in 1:Np) {
    for(c in 1:Nc) {
      sigma_condition[p][c] ~ cauchy(0, 1);
      mu_condition[p][c] ~ normal(prior_mean_condition, prior_sigma_condition);
    }
  }
  
  sigma_sample ~ cauchy(0, 1);
}


generated quantities {
  vector [Ng] Yhat_sample;
  
  for(g in 1:Ng) {
   Yhat_sample[g] = normal_rng(mu_sample[g], sigma_sample[g]);
  }
}
