data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ncol; // total number of columns in model matrix
  int<lower=1> NscaleBatches;
  matrix[Ntotal, Ncol] X; // model matrix
  real y[Ntotal]; // response variable t distributed
  int<lower=1, upper=NscaleBatches> NBatchMap[Ncol];
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression parameters
  real populationMean;//[NscaleBatches];
  real<lower=0.01> sigmaPop; // data standard deviation
  real<lower=0.01> sigmaRan[NscaleBatches]; // group level errors
  real<lower=1> nu; // normality parameter for t distribution or degree of freedom
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  // fitted value
  mu = X * betas;
  mu = mu + populationMean;
}
model {
  real sigmaRan_expanded[Ncol];
  //real populationMean_expanded[Ncol];
  nu ~ exponential(1/29.0);
  // using diffuse prior
  sigmaRan ~ cauchy(0, 2);
  sigmaPop ~ cauchy(0, 2);
  populationMean ~ cauchy(0, 1);
  // vector expansion by mapping to a larger vector/array
  sigmaRan_expanded = sigmaRan[NBatchMap];
  //populationMean_expanded = populationMean[NBatchMap];
  betas ~ normal(0, sigmaRan_expanded); //prior for the betas
  // likelihood function
  y ~ student_t(nu, mu, sigmaPop);
}
