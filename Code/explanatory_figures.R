library(ggplot2)
library(scales)
library(PNWColors)
theme_set(theme_classic())

e.axis <- c(-2,2,-1,1,0)
x.axis <- seq(0,1,length.out=1000)
dat <- data.frame(
  x=rep(x.axis, times=length(e.axis)),
  e=rep(e.axis, each=length(x.axis))
) %>%
  mutate(y=x^exp(-e))
ggplot(dat, aes(x=x,y=y,color=as.factor(e)))+
  geom_line(linewidth=1.5)+
  labs(x="Substrate percent cover",y=expression(paste(italic("O. suensonii")," percent occupancy")),color="Electivity")+
  scale_x_continuous(labels=percent)+
  scale_y_continuous(labels=percent)+
  scale_color_manual(values=pnw_palette("Sunset2",n=(N <- 7))[2:6])+
  theme(legend.position="top",text=element_text(size=22))

if (OUTPUT==TRUE)
  ggsave("fig/revision1/electivity_model.png", width=WID,height=WID, dpi=400)





if (FALSE) {
  ggplot(filter(dat,e%in%c(0)), aes(x=x,y=y,color=as.factor(e)))+
    geom_line(linewidth=1.5)+
    labs(x="Substrate percent cover",y=expression(paste(italic("O. suensonii")," percent occupancy")),color="Electivity")+
    scale_x_continuous(labels=percent)+
    scale_y_continuous(labels=percent)+
    scale_color_manual(values=pnw_palette("Sunset2",n=(N <- 7))[c(2:6)[3]])+
    theme(legend.position="top",text=element_text(size=22))
  ggsave("fig/revision1/electivity_model_1of3.png", width=WID,height=WID, dpi=400)
  ggplot(filter(dat,e%in%c(1,-1,0)), aes(x=x,y=y,color=as.factor(e)))+
    geom_line(linewidth=1.5)+
    labs(x="Substrate percent cover",y=expression(paste(italic("O. suensonii")," percent occupancy")),color="Electivity")+
    scale_x_continuous(labels=percent)+
    scale_y_continuous(labels=percent)+
    scale_color_manual(values=pnw_palette("Sunset2",n=(N <- 7))[c(2:6)[2:4]])+
    theme(legend.position="top",text=element_text(size=22))
  ggsave("fig/revision1/electivity_model_2of3.png", width=WID,height=WID, dpi=400)
  ggplot(dat, aes(x=x,y=y,color=as.factor(e)))+
    geom_line(linewidth=1.5)+
    labs(x="Substrate percent cover",y=expression(paste(italic("O. suensonii")," percent occupancy")),color="Electivity")+
    scale_x_continuous(labels=percent)+
    scale_y_continuous(labels=percent)+
    scale_color_manual(values=pnw_palette("Sunset2",n=(N <- 7))[c(2:6)])+
    theme(legend.position="top",text=element_text(size=22))
  ggsave("fig/revision1/electivity_model_3of3.png", width=WID,height=WID, dpi=400)
}


