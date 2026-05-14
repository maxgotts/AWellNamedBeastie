library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
library(tidyr)

theme_set(theme_classic()+theme(strip.background=element_blank(),text=element_text(family="Avenir")))

# Set this to TRUE if you want to overwrite the figure outputs
# Set this to FALSE if you don't want to save any of the figures
OUTPUT <- FALSE
if (FALSE) {
  OUTPUT <- TRUE
}

## SET UP ENVIRONMENT ##
# Set substrate names
ORDER <- c("Sponge","Octo", "Hydro", "Hexa")
ORDER_NM <- c("Sponge","Gorgonian", "Fire coral", "Scleractinian")

# Choose substrate colours
COLS <- c("#255f85","#c5283d", "#e9724c", "#481d24", "#848884")
names(COLS) <- ORDER_NM

# Set figure heights and style
WID <- 9.154167
HEI <- 5.820833

# Set working directory
setwd("~/Desktop/Research/BrittleStarResting/AWellNamedBeastie")

# Load data
star <- read.csv("Data/cleaned_brittle_stars.csv")

# Read model
fit <- readRDS("Model/electivity_output.rds")


## Make a function for generating posteriors
Posterior <- function(fit, pars, COLOURS, new_names, title="", credible=.83, total.alpha=0.9) {
  outputs <- as.data.frame(matrix(NA, nrow=0, ncol=3))
  colnames(outputs) <- c("x", "kde", "variable")
  draws <- fit |>
    as.data.frame(pars=pars) |>
    melt()
  for (par in pars) {
    drawsPar <- draws |>
      filter(variable==par) |>
      pull(value)
    xStart <- min(drawsPar)
    xEnd <- max(drawsPar)
    xPar <- seq(xStart, xEnd, length.out=500)
    limits83 <- quantile(drawsPar, probs=credible %>%
                           { (1-.)/2+c(0,1)*. }) |>
      as.vector()
    outputs <- drawsPar |>
      kde() |>
      predict(x=xPar) %>%
      {
        data.frame(x=xPar, kde=., variable=par)
      } |>
      mutate( # x >= limits83[1] & x <= limits83[2]
        credible = case_when(
          x >= limits83[1] & x <= limits83[2] ~ "Inner",
          .default = "Outer"
        )
      ) %>%
      {
        rbind(outputs, .)
      }
  }
  outputs$variable <- factor(new_names[outputs$variable,"new_names"], levels=new_names$new_names)
  
  plot <- outputs |>
    ggplot(aes(x=x,y=kde,fill=variable))+
    geom_area(color=NA, position="identity", alpha=0.5*total.alpha)+
    geom_area(color=NA, data=outputs |> filter(credible=="Inner"), color=NA, position="identity", alpha=0.5*total.alpha)+
    theme(
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.position="top",
      text=element_text(size=14)
    )+
    labs(x=title,y="",fill="")+
    scale_fill_manual(values=COLOURS)+
    scale_color_manual(values=COLOURS)
  return(plot)
}


## Make a look-up table to convert between the variable names and the pretty names
new_names <- data.frame(old_names=c("spongeElectivity","gorgonianElectivity","fireElectivity","scleractinianElectivity"), new_names=c("Sponges","Gorgonians","Fire corals","Scleractinians"))
rownames(new_names) <- new_names$old_names




## Visualize the posteriors
pars <- c("gorgonianElectivity","fireElectivity","spongeElectivity","scleractinianElectivity")
Posterior(fit, pars, as.character(COLS), title="Electivity", new_names=new_names)+
  geom_vline(xintercept=0, color="black", linetype="dashed")+
  scale_x_continuous(breaks=c(-1.5,-1,-.5,0,.5,1))+
  theme(strip.background=element_blank(),text=element_text(size=20))
if (OUTPUT==TRUE)
  ggsave("Figures/electivity_posteriors.png", width=WID,height=HEI, dpi=400)





## Create the cover-occupancy curves
set.seed(42) # Seed of 42 to ensure no one panics
ndraws <- 1000 # number of draws from the posterior
x.axis <- seq(0,1,length.out=100) # Cover axis goes from 0 to 1 in 1% increments
occ.abu <- as.data.frame(matrix(NA, nrow=0, ncol=4)) # This is going to hold the cover-occupancy curves
colnames(occ.abu) <- c("Substrate","Cover","Occupancy","Electivity")
posterior <- as.data.frame(fit, pars=pars)
for (i in 1:ndraws) {
  for (k in 1:length(ORDER)) {
    rows <- (nrow(occ.abu)+1):(nrow(occ.abu)+length(x.axis))
    occ.abu[rows,"Substrate"] <- ORDER_NM[k]
    occ.abu[rows,"Cover"] <- x.axis
    elect <- occ.abu[rows,"Electivity"] <- posterior[i,c("spongeElectivity","gorgonianElectivity","fireElectivity","scleractinianElectivity")[k]]
    occ.abu[rows,"Occupancy"] <- x.axis^(exp(-elect))
  }
}
occ.abu$Substrate <- factor(occ.abu$Substrate, levels=ORDER_NM)

# Create an estimate of the cover-abundance for each quadrat-substrate type using the real data
prop_by_cover <- star %>% 
  melt(variable.name="QuadSubstrate", value.name="X.cover", measure.vars=ORDER) %>% # Get the % cover of each substrate
  melt(variable.name="HostSubstrate",value.name="X.star",measure.vars=paste0("X.on.",tolower(ORDER))) %>% # Get the % occupancy of each substrate
  mutate(
    Host.is.Quad = HostSubstrate==paste0("X.on.",tolower(QuadSubstrate)),
    Substrate = QuadSubstrate
  ) %>% 
  filter(Host.is.Quad) %>% # Filter for only "% cover of [substrate]" and "% occupancy on [substrate]" for each substrate type
  select(SiteName, Substrate, X.cover, X.star)
prop_by_cover$Substrate <- factor(prop_by_cover$Substrate, levels=ORDER)
levels(prop_by_cover$Substrate) <- ORDER_NM

# Create figure
ggplot(occ.abu, aes(x=Cover,y=Occupancy,group=as.factor(Electivity),color=Substrate))+
  geom_line(alpha=0.1)+
  scale_x_continuous(labels=percent)+#,breaks=c(.25,.5,.75)
  scale_y_continuous(labels=percent)+#,breaks=c(.25,.5,.75)
  scale_color_manual(values=as.character(COLS))+
  scale_fill_manual(values=as.character(COLS))+
  geom_abline(slope=1, intercept=0, linetype="dashed",color="black")+#as.character(COLS[5])#data=data.frame(x=x.axis, y=x.axis), aes(x=x,y=y), inherit.aes=FALSE,
  theme(legend.position="none",text=element_text(size=30),strip.text=element_text(face="bold"),axis.text=element_text(size=20))+
  geom_point(data=prop_by_cover, aes(x=X.cover,y=X.star,fill=Substrate), color="white", pch=21, size=3, inherit.aes=FALSE)+
  facet_wrap(~Substrate)+
  theme(
    legend.position="none",
    strip.text=element_text(face="bold"),
    strip.background=element_blank(),
    axis.line = element_blank(),
    # axis.text.x=element_text(hjust=c(0,1))
    panel.spacing.x = unit(6.5, "mm"),
    plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt")
  )+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1)+
  labs(x="Substrate percent cover",y=expression(paste(italic("O. suensonii")," percent occupancy")))
if (OUTPUT==TRUE)
  ggsave("Figures/electivity_curves.png", width=WID,height=WID, dpi=400)





## Get 83% credible intervals
substrates <- c("spongeElectivity", "gorgonianElectivity", "fireElectivity", "scleractinianElectivity")
for (k in 1:4){
  substrate <- substrates[k]
  cat(substrate, "\n");
  print(round(quantile(as.data.frame(fit, pars=c(substrate))[,substrate], c(.5, .83 %>% { (1-.)/2+c(0,1)*. }) %>% sort),3));
  if (k!=4) cat("\n")
}

## Get P(par>0)
as.data.frame(fit, pars=substrates) %>% 
  melt %>% 
  mutate(above.0=value>0) %>%
  group_by(above.0,variable) %>%
  tally %>%
  mutate(above.0=c("Below","Above")[as.numeric(above.0)+1]) %>% 
  pivot_wider(id_cols=variable,names_from=above.0,values_from=n,values_fill=0) %>% 
  mutate(Total=Below+Above,X.below=100*Below/Total,X.above=100*Above/Total)




## Get % individuals above or below line
prop_by_cover %>% 
  mutate(
    Delta = X.star - X.cover,
    WhereOnLine = ifelse(Delta > 0,"Above","Below")
  ) %>%
  group_by(Substrate, WhereOnLine) %>% 
  tally() %>% 
  group_by(Substrate) %>% 
  mutate(Proportion=100*n/sum(n),Total=sum(n)) %>% 
  ungroup %>% 
  rename(Subset=n) %>% 
  select(Substrate,WhereOnLine,Subset,Total,Proportion)
