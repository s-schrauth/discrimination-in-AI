//
// Adapted code from https://github.com/mkusner/counterfactual-fairness
//

// Input data
data {
  
  int<lower = 0> N;          // number of cases
  int<lower = 0> K;          // number of protected variables
  matrix[N, K]   prot;       // protected variables
  vector[N]      age; 
  array[N] int   priors; 
  array[N] int   crime; 
  array[N] int   juvmisd; 
  array[N] int   juvfel;
  array[N] int   compas;
  
}

// Help vectors for modelling process
transformed data {
  
 vector[K] zero_K;
 vector[K] one_K;
 
 zero_K = rep_vector(0,K);
 one_K = rep_vector(1,K);

}

// Parameter definition
parameters {
  
  // latent variables
  vector[N] u_d;
  vector[N] u_j;
  
  // weights of observable variables
  real priors_bias;
  real priors_ud;
  real priors_age;
  real crime_bias;
  real crime_ud;
  real crime_age;
  real juvmisd_bias;
  real juvmisd_uj;
  real juvmisd_age;
  real juvfel_bias;
  real juvfel_uj;
  real juvfel_age;
  real compas_bias;
  real compas_ud;
  real compas_age;
  
  // weight vectors of protected variable set
  vector[K] priors_prot;
  vector[K] crime_prot;
  vector[K] juvmisd_prot;
  vector[K] juvfel_prot;
  vector[K] compas_prot;
  
}

// Prior and model definition
model {
  
  // Priors
  u_d           ~ normal(0,1);
  u_j           ~ normal(0,1);
  
  priors_bias   ~ normal(0,1);
  priors_ud     ~ normal(0,1);
  priors_age    ~ normal(0,1);
  crime_bias    ~ normal(0,1);
  crime_ud      ~ normal(0,1);
  crime_age     ~ normal(0,1);
  juvmisd_bias  ~ normal(0,1);
  juvmisd_uj    ~ normal(0,1);
  juvmisd_age   ~ normal(0,1);
  juvfel_bias   ~ normal(0,1);
  juvfel_uj     ~ normal(0,1);
  juvfel_age    ~ normal(0,1);
  compas_bias   ~ normal(0,1);
  compas_ud     ~ normal(0,1);
  compas_age    ~ normal(0,1);
  
  priors_prot   ~ normal(zero_K, one_K);
  crime_prot    ~ normal(zero_K, one_K);
  juvmisd_prot  ~ normal(zero_K, one_K);
  juvfel_prot   ~ normal(zero_K, one_K);
  compas_prot   ~ normal(zero_K, one_K);
  
  // Model
  priors  ~ poisson(exp(priors_bias + priors_ud * u_d + priors_age * age + prot * priors_prot));
  crime   ~ bernoulli(Phi(crime_bias + crime_ud * u_d + crime_age * age + prot * crime_prot));
  juvmisd ~ poisson(exp(juvmisd_bias + juvmisd_uj * u_j + juvmisd_age * age + prot * juvmisd_prot));
  juvfel  ~ poisson(exp(juvfel_bias + juvfel_uj * u_j + juvfel_age * age + prot * juvfel_prot));
  
  compas  ~ bernoulli_logit(compas_bias + compas_ud * u_d + compas_age * age + prot * compas_prot);
}
