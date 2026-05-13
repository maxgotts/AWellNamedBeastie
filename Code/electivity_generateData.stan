data {
  int nQuad; // Number of quadrats
  vector[nQuad] gorgonianAbundance; // Quadrat cover of gorgonians
  vector[nQuad] fireAbundance; // Quadrat cover of fire corals
  vector[nQuad] spongeAbundance; // Quadrat cover of sponges
  vector[nQuad] scleractinianAbundance; // Quadrat cover of scleractinians
  array[nQuad] int totalSuensonii; // Total number of suensonii in quadrat
  real gorgonianElectivity; // Quadrat cover of gorgonians
  real fireElectivity; // Quadrat cover of fire corals
  real spongeElectivity; // Quadrat cover of sponges
  real scleractinianElectivity; // Quadrat cover of scleractinians
}

generated quantities {
  real positiveGorgonianElectivity = exp(-gorgonianElectivity); // Rescaled gorgonian electivity
  real positiveFireElectivity = exp(-fireElectivity); // Rescaled fire coral electivity
  real positiveSpongeElectivity = exp(-spongeElectivity); // Rescaled sponge electivity
  real positiveScleractinianElectivity = exp(-scleractinianElectivity); // Rescaled scleractinian electivity
  
  array[nQuad] int gorgonianSuensonii = binomial_rng(totalSuensonii, gorgonianAbundance^positiveGorgonianElectivity);  // Total number of suensonii on gorgonians
  array[nQuad] int fireSuensonii = binomial_rng(totalSuensonii, fireAbundance^positiveFireElectivity); // Total number of suensonii on fire corals
  array[nQuad] int spongeSuensonii = binomial_rng(totalSuensonii, spongeAbundance^positiveSpongeElectivity); // Total number of suensonii on sponges
  array[nQuad] int scleractinianSuensonii = binomial_rng(totalSuensonii, scleractinianAbundance^positiveScleractinianElectivity); // Total number of suensonii on scleractinians
}




