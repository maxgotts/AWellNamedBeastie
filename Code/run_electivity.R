library(rstan)
library(bayesplot)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(reshape2)

# Set the theme for this
theme_set(theme_classic()+theme(text=element_text(family="Avenir")))

# Initialize Stan with parallelization 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set this to TRUE if you want to overwrite the model outputs
# Set this to FALSE if you don't want to write out any of the results
OUTPUT <- FALSE
if (FALSE) {
  OUTPUT <- TRUE
}


## INITIALIZE DATASETS ##
# Load dataset
setwd("~/Desktop/Research/BrittleStarResting/AWellNamedBeastie")
star <- read.csv("Data/brittle_stars.csv")

# Turn data into proportions
star$TotalCover <- star$Sponge.. + star$Hydro..+ star$Octo.. + star$Hexa..
star$Sponge <- star$Sponge../star$TotalCover
star$Hydro <- star$Hydro../star$TotalCover
star$Octo <- star$Octo../star$TotalCover
star$Hexa <- star$Hexa../star$TotalCover
star$Sponge.. <- NULL
star$Hydro.. <- NULL
star$Octo.. <- NULL
star$Hexa.. <- NULL
star$Total.star <- star$Brittle.on.sponge+star$Brittle.on.hydro+star$Brittle.on.octo+star$Brittle.on.hexa
star$X.on.sponge <- star$Brittle.on.sponge/star$Total.star
star$X.on.hydro <- star$Brittle.on.hydro/star$Total.star
star$X.on.octo <- star$Brittle.on.octo/star$Total.star
star$X.on.hexa <- star$Brittle.on.hexa/star$Total.star

write.csv(star,"Data/cleaned_brittle_stars.csv",row.names=FALSE)


model <- stan_model("Code/electivity_model.stan")

# Assemble real data
real_data <- list(
  nQuad = nrow(star),
  gorgonianAbundance = star$Octo,
  fireAbundance = star$Hydro,
  spongeAbundance = star$Sponge,
  scleractinianAbundance = star$Hexa,
  totalSuensonii = star$Total.star,
  gorgonianSuensonii = star$Brittle.on.octo,
  fireSuensonii = star$Brittle.on.hydro,
  spongeSuensonii = star$Brittle.on.sponge,
  scleractinianSuensonii = star$Brittle.on.hexa
)

fit <- sampling(
  object = model,
  data = real_data,
  chains = 4,
  seed = 42,
  iter = (itr <- 12000),
  warmup = (warm <- itr/2),
  control = options(max_treedepth=12, adapt_delta=0.9999)
) # You may have to run this twice, depending on your Stan set-up

# Save (or read in) model output
model_path <- paste0("Model/electivity_output.rds")
if (FALSE)
  saveRDS(fit, model_path)
if (FALSE)
  fit <- readRDS(model_path)



## CHECK FIT QUALITY
pars <- c("gorgonianElectivity","fireElectivity","spongeElectivity","scleractinianElectivity")

fit %>% mcmc_trace(pars=pars)
fit %>% mcmc_dens_overlay(pars=pars)
fit %>% mcmc_acf(pars=pars)
fit %>% pairs(pars=pars)

rhat(fit) %T>%
  print %>%
  range %>%
  round(5) %>% 
  paste0(collapse=" - ") %>%
  paste0("\n\nRhat range = ",.) %>%
  cat




