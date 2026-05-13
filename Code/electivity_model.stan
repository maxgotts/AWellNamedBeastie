data {
  int nQuad; // Number of quadrats
  vector[nQuad] gorgonianAbundance; // Quadrat cover of gorgonians
  vector[nQuad] fireAbundance; // Quadrat cover of fire corals
  vector[nQuad] spongeAbundance; // Quadrat cover of sponges
  vector[nQuad] scleractinianAbundance; // Quadrat cover of scleractinians
  array[nQuad] int totalSuensonii; // Total number of suensonii in quadrat
  array[nQuad] int gorgonianSuensonii; // Total number of suensonii on gorgonians
  array[nQuad] int fireSuensonii; // Total number of suensonii on fire corals
  array[nQuad] int spongeSuensonii; // Total number of suensonii on sponges
  array[nQuad] int scleractinianSuensonii; // Total number of suensonii on scleractinians
}

transformed data {
  real epsilon = 1e-6;
  vector[nQuad] cushionGorgonianAbundance = fmin(1-epsilon, gorgonianAbundance+epsilon);
  vector[nQuad] cushionFireAbundance = fmin(1-epsilon,fireAbundance+epsilon);
  vector[nQuad] cushionSpongeAbundance = fmin(1-epsilon, spongeAbundance+epsilon);
  vector[nQuad] cushionScleractinianAbundance = fmin(1-epsilon, scleractinianAbundance+epsilon);
}

parameters {
  real<lower=0> unscaledGorgonianElectivity; // Unscaled gorgonian electivity 
  real<lower=0> unscaledFireElectivity; // Unscaled fire coral electivity
  real<lower=0> unscaledSpongeElectivity; // Unscaled sponge electivity
  real<lower=0> unscaledScleractinianElectivity; // Unscaled scleractinian electivity
}
  
model {
  log(unscaledGorgonianElectivity) ~ normal(0,2);
  log(unscaledFireElectivity) ~ normal(0,2);
  log(unscaledSpongeElectivity) ~ normal(0,2);
  log(unscaledScleractinianElectivity) ~ normal(0,2);
  // unscaledGorgonianElectivity ~ lognormal(0,2);
  // unscaledFireElectivity ~ lognormal(0,2);
  // unscaledSpongeElectivity ~ lognormal(0,2);
  // unscaledScleractinianElectivity ~ lognormal(0,2);
  
  gorgonianSuensonii ~ binomial(totalSuensonii, cushionGorgonianAbundance^unscaledGorgonianElectivity);
  fireSuensonii ~ binomial(totalSuensonii, cushionFireAbundance^unscaledFireElectivity);
  spongeSuensonii ~ binomial(totalSuensonii, cushionSpongeAbundance^unscaledSpongeElectivity);
  scleractinianSuensonii ~ binomial(totalSuensonii, cushionScleractinianAbundance^unscaledScleractinianElectivity);
}

generated quantities {
  real gorgonianElectivity = -log(unscaledGorgonianElectivity); // Rescaled gorgonian electivity
  real fireElectivity = -log(unscaledFireElectivity); // Rescaled fire coral electivity
  real spongeElectivity = -log(unscaledSpongeElectivity); // Rescaled sponge electivity
  real scleractinianElectivity = -log(unscaledScleractinianElectivity); // Rescaled scleractinian electivity
}


