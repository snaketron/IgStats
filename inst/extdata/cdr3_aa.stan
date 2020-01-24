data {
  int <lower=0> Na;                 // number of amino acids => 20
  int <lower=0> Ng;                 // number of samples
  int <lower=0> Nc;                 // number of conditions
  int <lower=0> Y [Ng, Na];         // CDR3 length
  int <lower=0> C [Ng];             // condition ID
}


parameters {
  simplex [Na] mu_sample [Ng];
  vector <lower=0> [Na] alpha [Nc];
}


model {
  int Yvec [Na];
  for(s in 1:Ng) {
    Y[s, ] ~ multinomial(mu_sample[s]);
    mu_sample[s] ~ dirichlet(alpha[C[s]]);
  }
  for(j in 1:Nc) {
    alpha[j] ~ normal(1, 10);
  }
}

generated quantities {
  int <lower=0> Yhat_sample [Ng, Na];
  real <lower=0> Yhat_sample_real [Ng, Na];
  vector [Na] mu_condition [Nc];

  for(s in 1:Ng) {
    Yhat_sample[s, ] = multinomial_rng(mu_sample[s], sum(Y[s,]));
    Yhat_sample_real[s, ] = Yhat_sample[s, ];
    for(a in 1:Na) {
      Yhat_sample_real[s, a] = Yhat_sample_real[s, a] / sum(Y[s,]);
    }
  }

  for(c in 1:Nc) {
    // softmax
    mu_condition[c] = exp(log(alpha[c]))/sum(exp(log(alpha[c])));
  }
}
