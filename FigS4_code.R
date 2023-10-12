library(tidyverse)
source("./functions.R")

fniche<-tibble("spi"=stnames,"fundniche"=fundniche)
int<-read_csv("data4_interactions.csv") %>% left_join(fniche,by="spi") %>% 
  mutate(Pratio=(10^Pi)/fundniche)

g<-ggplot(data=int,aes(x=10^Pi/1000000,y=Cij))+
  geom_point(pch=21,size=3,alpha=0.8,aes(fill=spj))+
  theme_bw(base_size=18)+
  scale_x_continuous(limits=c(0,NA),expand=c(0,0.1),labels=scalefun)+
  scale_y_continuous(limits=c(-0.22,0.13),expand=c(0,0))+
  stat_smooth(method="lm",color="black",size=1.2)+
  scale_fill_manual(name="Affecting species",values=stcolors)+
  geom_text(data=distinct(int,spi),aes(x=Inf,y=0.1,label=spi,hjust=1.2),size=6)+
  theme(
    panel.border=element_blank(),
    axis.line=element_line(color="black"),
    aspect.ratio=0.7,
    panel.grid=element_blank(),
    axis.title=element_blank(),
    axis.text=element_text(color="black"),
    strip.text=element_blank(),
    strip.background=element_blank(),
    legend.position=c(0.9,0.2)
  )+
  facet_wrap(~spi,ncol=4,scales="free")

win.graph(18,8)
g 

ggsave("FigS4.pdf")
