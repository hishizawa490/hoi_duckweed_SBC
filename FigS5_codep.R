library(tidyverse)
library(agricolae)
library(PTXQC)
library(patchwork)
setwd("C:/Users/uadgw/デスクトップ/WorkingFiles/モデル群集/syncomR/deposit")
source("./functions.R")
int<-read_csv("data4_interactions.csv")
pred3<-read_csv("data6_abundance_prediction_trio.csv")
predint3<-extract_int(pred3) %>% mutate(type=3)
pred2<-read_csv("data7_abundance_prediction_pair.csv")
predint2<-extract_int(pred2) %>% mutate(type=2)
int[,5:7]<-10^int[5:7]/1000000
predint3[,6:8]<-10^predint3[6:8]/1000000
predint2[,6:8]<-10^predint2[6:8]/1000000

max<-int %>% group_by(spi) %>% summarize(maxobs=max(Pi.j))
g1<-ggplot(data=int)+
  geom_segment(aes(x=nspecies-1,xend=nspecies,y=Pi,yend=Pi.j,color=str_c("+ ",spj)),size=0.55)+
  geom_point(aes(x=nspecies,y=Pi.j),size=1.6,shape=21,fill="white",stroke=0.8)+
  geom_point(data=filter(int,nspecies==2),aes(x=1,y=Pi),size=1.6,shape=21,fill="white",stroke=0.8)+
  geom_point(data=max,aes(x=1,y=maxobs*1.08),fill="transparent",color="transparent")+
  scale_color_manual(values=stcolors,name="")+
  scale_x_continuous(breaks=seq(1,max(int$nspecies),1))+
  scale_y_continuous(labels=scalefun,limits=c(0,NA),breaks=breaks_fun2,expand=c(0,0))+
  theme_bw(base_size=16)+
  guides(color=guide_legend(override.aes=list(size=1.4)))+
  theme(
    legend.position="none",
    axis.line=element_line(color="black"),
    axis.title=element_blank(),
    axis.text.x=element_text(size=14,color="black",margin=margin(0.1,0.1,0.1,0.1,"lines")),
    axis.text.y=element_text(size=14,color="black",margin=margin(0.1,0.1,0.1,0.1,"lines")),
    axis.ticks.length=unit(0.2,"lines"),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    aspect.ratio=0.85,
    strip.text=element_text(size=16,margin=margin(0,0,0.3,0,"lines")),
    strip.background=element_blank(),
    panel.spacing.x=unit(0.3,"lines"),
    panel.spacing.y=unit(0.3,"lines")
  )+
  coord_cartesian(ylim=c(-0.015,NA))+
  facet_wrap(~spi,nrow=2,scales="free")

max<-predint3 %>% group_by(spi) %>% summarize(maxobs=max(Pi.j))
g2<-ggplot(data=predint3)+
  geom_segment(aes(x=nspecies-1,xend=nspecies,y=Pi,yend=Pi.j,color=str_c("+ ",spj)),size=0.55)+
  geom_point(aes(x=nspecies,y=Pi.j),size=1.6,shape=21,fill="white",stroke=0.8)+
  geom_point(data=filter(predint3,nspecies==2),aes(x=1,y=Pi),size=1.6,shape=21,fill="white",stroke=0.8)+
  geom_point(data=max,aes(x=1,y=maxobs*1.08),fill="transparent",color="transparent")+
  scale_color_manual(values=stcolors,name="")+
  scale_x_continuous(breaks=seq(1,max(int$nspecies),1))+
  scale_y_continuous(labels=scalefun,limits=c(0,NA),breaks=breaks_fun2,expand=c(0,0))+
  theme_bw(base_size=16)+
  guides(color=guide_legend(override.aes=list(size=1.4)))+
  theme(
    legend.position="none",
    axis.line=element_line(color="black"),
    axis.title=element_blank(),
    axis.text.x=element_text(size=14,color="black",margin=margin(0.1,0.1,0.1,0.1,"lines")),
    axis.text.y=element_text(size=14,color="black",margin=margin(0.1,0.1,0.1,0.1,"lines")),
    axis.ticks.length=unit(0.2,"lines"),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    aspect.ratio=0.85,
    strip.text=element_text(size=16,margin=margin(0,0,0.3,0,"lines")),
    strip.background=element_blank(),
    panel.spacing.x=unit(0.3,"lines"),
    panel.spacing.y=unit(0.3,"lines")
  )+
  coord_cartesian(ylim=c(-0.015,NA))+
  facet_wrap(~spi,nrow=2,scales="free")

max<-predint3 %>% group_by(spi) %>% summarize(maxobs=max(Pi.j))
g3<-ggplot(data=predint2)+
  geom_segment(aes(x=nspecies-1,xend=nspecies,y=Pi,yend=Pi.j,color=str_c("+ ",spj)),size=0.55)+
  geom_point(aes(x=nspecies,y=Pi.j),size=1.6,shape=21,fill="white",stroke=0.8)+
  geom_point(data=filter(predint2,nspecies==2),aes(x=1,y=Pi),size=1.6,shape=21,fill="white",stroke=0.8)+
  geom_point(data=max,aes(x=1,y=maxobs*1.08),fill="transparent",color="transparent")+
  scale_color_manual(values=stcolors,name="")+
  scale_x_continuous(breaks=seq(1,max(int$nspecies),1))+
  scale_y_continuous(labels=scalefun,limits=c(0,NA),breaks=breaks_fun2,expand=c(0,0.1))+
  theme_bw(base_size=16)+
  guides(color=guide_legend(override.aes=list(size=1.4)))+
  theme(
    legend.position="none",
    axis.line=element_line(color="black"),
    axis.title=element_blank(),
    axis.text.x=element_text(size=14,color="black",margin=margin(0.1,0.1,0.1,0.1,"lines")),
    axis.text.y=element_text(size=14,color="black",margin=margin(0.1,0.1,0.1,0.1,"lines")),
    axis.ticks.length=unit(0.2,"lines"),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    aspect.ratio=0.85,
    strip.text=element_text(size=16,margin=margin(0,0,0.3,0,"lines")),
    strip.background=element_blank(),
    panel.spacing.x=unit(0.3,"lines"),
    panel.spacing.y=unit(0.3,"lines")
  )+
  coord_cartesian(ylim=c(-0.015,NA))+
  facet_wrap(~spi,nrow=2,scales="free")

win.graph(30,40)
g1/g2/g3



ggsave("FigS5.pdf")   
dev.off()
