library(rstan)
library(bayesplot)
library(tidybayes)
library(ggridges)
library(gridExtra)
library(ks)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(broom)
library(brms)

# Initialize stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Initialize ggplot2
theme_set(theme_classic()+theme(text=element_text(size=14),strip.background=element_blank())+theme(text=element_text(family="Avenir")))



# Set working directory
setwd("~/Desktop/Research/BrittleStarResting/Code/ELECTIVITY")


### GENERATE SIMULATION ###
# Load model
simulation_model <- stan_model(file = "electivity_generateData.stan")
analysis_model <- stan_model(file = "electivity_model.stan")

# Parameter values
set.seed(42)
retrodiction <- as.data.frame(matrix(NA,nrow=0,ncol=9))
NUMBER_OF_RUNS <- 50
for (modelrun in 1:NUMBER_OF_RUNS) {
  simulation_parameters <- list(
    nQuad=(nQuad<-32), # Number of quadrats
    gorgonianAbundance=runif(nQuad, min=0, max=1), # Quadrat cover of gorgonians
    fireAbundance=runif(nQuad, min=0, max=1), # Quadrat cover of fire corals
    spongeAbundance=runif(nQuad, min=0, max=1), # Quadrat cover of sponges
    scleractinianAbundance=runif(nQuad, min=0, max=1), # Quadrat cover of scleractinians -- intentionally empty!
    totalSuensonii=rpois(nQuad, lambda=9), # Total number of suensonii in quadrat
    gorgonianElectivity=rnorm(1,0,1), # Quadrat cover of gorgonians
    fireElectivity=rnorm(1,0,1), # Quadrat cover of fire corals
    spongeElectivity=rnorm(1,0,1), # Quadrat cover of sponges
    scleractinianElectivity=rnorm(1,0,1) # Quadrat cover of scleractinians
  )
  sum_to_one <- simulation_parameters$gorgonianAbundance+simulation_parameters$fireAbundance+simulation_parameters$spongeAbundance+simulation_parameters$scleractinianAbundance
  simulation_parameters$gorgonianAbundance <- simulation_parameters$gorgonianAbundance/sum_to_one
  simulation_parameters$fireAbundance <- simulation_parameters$fireAbundance/sum_to_one
  simulation_parameters$spongeAbundance <- simulation_parameters$spongeAbundance/sum_to_one
  simulation_parameters$scleractinianAbundance <- simulation_parameters$scleractinianAbundance/sum_to_one
  
  # Run model
  sim <- sampling(
    object = simulation_model,
    data = simulation_parameters,
    chains = 1,
    seed = 42,
    iter = 1,
    algorithm = "Fixed_param"
  )
  
  # Save (or read in) simulation model output
  simulation_model_path <- paste0("Retrodiction_output/simulation_run_",modelrun,".rds")
  # if (FALSE)
    saveRDS(sim, simulation_model_path)
  if (FALSE)
    sim <- readRDS(simulation_model_path)
  
  # Save (or read in) simulated data output
  simulation_dat <- sim %>% 
    as.data.frame(pars=c("gorgonianSuensonii","fireSuensonii","spongeSuensonii","scleractinianSuensonii")) %>% 
    reform() %>% 
    mutate(index=NULL)
  
  # simulation_dataframe_path <- paste0("Retrodiction_output/electivity_retrodiction_d",modelrun,".csv") 
  # if (FALSE)
  #   write.csv(simulation_dat, simulation_dataframe_path, row.names=F)
  # if (FALSE)
  #   simulation_dat <- read.csv(simulation_dataframe_path)
  
  ### ANALYSE SIMULATION ###
  # Assemble data
  simulated_data <- list(
    nQuad = nQuad,
    gorgonianAbundance = simulation_parameters$gorgonianAbundance,
    fireAbundance = simulation_parameters$fireAbundance,
    spongeAbundance = simulation_parameters$spongeAbundance,
    scleractinianAbundance = simulation_parameters$scleractinianAbundance,
    totalSuensonii = simulation_dat$gorgonianSuensonii+simulation_dat$fireSuensonii+simulation_dat$spongeSuensonii+simulation_dat$scleractinianSuensonii,
    gorgonianSuensonii = simulation_dat$gorgonianSuensonii,
    fireSuensonii = simulation_dat$fireSuensonii,
    spongeSuensonii = simulation_dat$spongeSuensonii,
    scleractinianSuensonii = simulation_dat$scleractinianSuensonii
  )
  
  # Run model
  fit <- sampling(
    object = analysis_model,
    data = simulated_data,
    chains = 4,
    seed = 42,
    iter = (itr <- 12000),
    warmup = (warm <- itr/2),
    control = options(max_treedepth=12, adapt_delta=0.9999)
  )
  
  # Save (or read in) model output
  # model_path <- paste0("Model_outputs/model1_m1v",modelversion,"_on_model1_data1_r",modelrun,".rds") 
  analysis_model_path <- paste0("Retrodiction_output/analysis_run_",modelrun,".rds")
  # if (FALSE)
    saveRDS(fit, analysis_model_path)
  if (FALSE)
    fit <- readRDS(model_path)
  
  # Name parameters of interest
  pars <- c("gorgonianElectivity","fireElectivity","spongeElectivity","scleractinianElectivity")
  
  # Test model fit
  if (FALSE) {
    fit %>% mcmc_trace(pars=pars)
    fit %>% mcmc_dens_overlay(pars=pars)
    fit %>% mcmc_acf(pars=pars)
    fit %>% pairs(pars=pars)
  }
  
  ### EVALUATE BIAS ###
  # Visualize bias
  new_results <- fit %>% 
    as.data.frame(pars=pars) %>% 
    melt() %>% 
    group_by(variable) %>% 
    summarize(
      lower = get83(value)[1],
      upper = get83(value)[2],
      med = median(value)
    ) %>%
    {
      cbind(., real=c(
        simulation_parameters$gorgonianElectivity,
        simulation_parameters$fireElectivity,
        simulation_parameters$spongeElectivity,
        simulation_parameters$scleractinianElectivity
      ))
    } %>% 
    mutate(
      real.in.83 = real <= upper & real >= lower,
      real.minus.med = real-med,
      percent.error = 100*(med-real)/real,
      run_id = modelrun
    )
  retrodiction <- rbind(retrodiction, new_results)

  # # Get data out
  # posterior_dat <- fit %>%
  #   as.data.frame(pars=pars) %>% 
  #   melt()
  # simulation_dat <- data.frame(
  #   gorgonianElectivity = simulation_parameters$gorgonianElectivity,
  #   fireElectivity = simulation_parameters$fireElectivity,
  #   spongeElectivity = simulation_parameters$spongeElectivity,
  #   scleractinianElectivity = simulation_parameters$scleractinianElectivity
  # ) %>% 
  #   melt()
  # 
  # # Visualize results
  # ggplot()+
  #   geom_density(data=posterior_dat, aes(x=value,fill=variable), color=NA, alpha=0.7)+
  #   geom_vline(data=simulation_dat, aes(xintercept=value, color=variable))+
  #   facet_wrap(~variable)+
  #   theme(legend.position="none")
  # # Posterior(fit, pars=pars, new_names = ORDER_NM)+
  # #   geom_vline(data=simulation_dat, aes(xintercept=value, color=variable))+
  # #   facet_wrap(~variable)+
  # #   theme(legend.position="none")
}

retrodiction %>% 
  group_by(variable, real.in.83) %>% 
  tally %>% 
  # filter(!real.in.83) %>%
  mutate(percent.wrong=100*n/NUMBER_OF_RUNS) %>% View
# if (FALSE)
  write.csv(retrodiction,"Retrodiction_output/retrodiction_results.csv",row.names=F)
if (FALSE)
  retrodiction <- read.csv("Retrodiction_output_RUN1_50times/retrodiction_results.csv")

retrodiction_annot <- retrodiction %>% 
  mutate(
    real_above_0 = lower > 0 & upper > 0,
    real_below_0 = lower < 0 & upper < 0,
    real_at_0 = lower < 0 & upper > 0,
    posterior = case_when(
      real_above_0~"positive",
      real_below_0~"negative",
      real_at_0~"zero",
      .default=NA
    ),
    input = ifelse(real>0,"positive","negative")
  ) %>% {
    .$variable <- factor(.$variable,levels=paste0(c("sponge","gorgonian","fire","scleractinian"),"Electivity"))
    levels(.$variable) <- c("Sponge","Gorgonian","Fire coral","Scleractinian")
    .
  }

retrodiction_annot %>% 
  ggplot(aes(x=real,y=med,color=posterior))+
  # facet_wrap(~variable)+
  geom_point()+
  geom_linerange(aes(ymin=lower,ymax=upper))+
  geom_abline(intercept=0,slope=1,linetype="dashed",color="#666666")+
  geom_vline(xintercept=0,linetype="dashed",color="#666666")+
  geom_hline(yintercept=0,linetype="dashed",color="#666666")+
  # labs(x="Real parameter value",y="Estimated parameter value (with 83% CI)")+
  labs(x="Real electivity",y="Estimated electivity")+
  scale_colour_manual(values=c("#FF6542","#648E90","#272727"))+
  theme(legend.position="none")+
  geom_smooth(aes(x=real,y=med),method="lm",inherit.aes = F,color="#272727",fill="#272727",linetype="solid")
if (FALSE)
  ggsave("Retrodiction_output_RUN1_50times/retrodiction_output2.jpeg",width=6,height=6,units="in")
  
  
retrodiction_annot %>% 
  # group_by(variable) %>% 
  do(tidy(lm(med~real,.))) %>% 
  filter(term=="real") %>% 
  select(
    # variable,
    estimate,
    std.error
  ) %>% 
  rename(
    # Substrate=variable,
    Slope=estimate,
    SE=std.error
  ) %>% 
  mutate(
    Lower=Slope-SE,
    Upper=Slope+SE
  )

lm.fit <- retrodiction_annot %>% 
  brm(med~real,data=.,iter=2000,warmup=1000,seed=42)
lm.fit %>% as.data.frame %>% View
retrodiction_annot %>%
  mutate(diff=real-med) %>% 
  # ggdens(diff)
  pull(diff) %>% sd

retrodiction_annot %>%
  # filter(posterior!="zero") %>% 
  group_by(
    # variable,
    input,
    posterior
  ) %>% 
  tally %>% 
  group_by(
    # variable,
    input
  ) %>% 
  mutate(
    RECALL_of_input = 100*n/sum(n) # TP out of total T
  ) %>% 
  ungroup %>% 
  filter(input==posterior) %>% 
  select(
    # variable,
    posterior,
    RECALL_of_input
  ) %>% 
  rename(
    input=posterior,
    # Substrate=variable
  ) %>%
  pivot_wider(
    names_from=input,
    values_from=RECALL_of_input
  ) %>% 
  rename(
    RECALL_of_positive=positive,
    RECALL_of_negative=negative
  )

retrodiction_annot %>% 
  # filter(posterior!="zero") %>%
  group_by(
    # variable,
    input,
    posterior
  ) %>% 
  tally %>% 
  group_by(
    # variable,
    posterior
  ) %>% 
  mutate(
    PRECISION_of_input = 100*n/sum(n) # TP out of total P
  ) %>% 
  ungroup %>% 
  filter(input==posterior) %>% 
  select(
    # variable,
    posterior,
    PRECISION_of_input
  ) %>% 
  rename(
    input=posterior,
    # Substrate=variable
  ) %>%
  pivot_wider(
    names_from=input,
    values_from=PRECISION_of_input
  ) %>% 
  rename(
    PRECISION_of_positive=positive,
    PRECISION_of_negative=negative
  )


retrodiction_annot %>% 
  # filter(variable=="Sponge") %>%
  group_by(
    # variable,
    input,
    posterior
  ) %>% 
  tally %>% 
  # group_by(variable) %>%
  mutate(x="brittle star") %>% group_by(x) %>% 
  mutate(total=sum(n)) %>% 
  filter(input==posterior) %>% 
  reframe(
    ACCURACY=100*sum(n)/total # T out of total
  ) %>% 
  # rename(Substrate=variable) %>% 
  distinct
