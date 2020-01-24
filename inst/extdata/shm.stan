data {
  int<lower=0> Ny;      // number of sequences (observations)
  int<lower=0> Ns;      // number of subjects
  int<lower=0> Nc;      // number of conditions
  int<lower=0> Nr;      // number of regions
  int<lower=0> N [Ny];  // sequence length
  int<lower=0> Y [Ny];  // number of mutations in sequence
  int S [Ny];           // sample IDs
  int R [Ny];           // sequence region ID
  int C [Ns];           // condition IDs
}


parameters {
  vector <lower=0, upper=1> [Nr] mu_subject [Ns];
  vector <lower=0, upper=1> [Nr] mu_condition [Nc];
  vector <lower=0> [Nc] phi;
  vector <lower=0> [Nc] tau;
}


model {
  
  for(i in 1:Ny) {
    Y[i] ~ binomial(N[i], mu_subject[S[i]][R[i]]);
  }
  
  for(r in 1:Nr) {
    for(s in 1:Ns) {
      mu_subject[s][r] ~ beta(mu_condition[C[s]][r]*phi[C[s]], phi[C[s]] - mu_condition[C[s]][r]*phi[C[s]]);
    }
  }
  
  for(c in 1:Nc) {
    // pareto II
    phi[c] ~ exponential(tau[c]);
    tau[c] ~ gamma(3, 0.1);
    for(r in 1:Nr) {
      mu_condition[c][r] ~ beta(1, 1);
    }
  }
}
