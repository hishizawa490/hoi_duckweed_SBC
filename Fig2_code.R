library(tidyverse)
source("./functions.R")
int<-read_csv("data4_interactions.csv")
int[,5:7]<-10^int[5:7]/1000000

max<-int %>% group_by(spi) %>% summarize(maxobs=max(Pi.j))

g<-ggplot(data=int)+
  geom_segment(aes(x=nspecies-1,xend=nspecies,y=Pi,yend=Pi.j,color=str_c(spj)),size=0.55)+
  geom_point(aes(x=nspecies,y=Pi.j),size=1.4,shape=21,fill="white",stroke=0.6)+
  geom_point(data=filter(int,nspecies==2),aes(x=1,y=Pi),size=1.4,shape=21,fill="white",stroke=0.8)+
  geom_point(data=max,aes(x=1,y=maxobs*1.08),fill="transparent",color="transparent")+
  scale_color_manual(values=stcolors,name="")+
  scale_x_continuous(breaks=seq(1,max(int$nspecies),1),name="Number of species")+
  scale_y_continuous(labels=scalefun,limits=c(0,NA),breaks=breaks_fun2,expand=c(0,0),
                     name=expression(paste(" Abundance (×",{10^6}," cfu ",{mg^-1},"plant)")))+
  theme_bw(base_size=20)+
  guides(color=guide_legend(override.aes=list(size=1.4)))+
  theme(
    legend.position=c(0.92,0.5),
    legend.justification=c(1,1),
    legend.background=element_rect(fill=NA,color=NA),
    legend.text=element_text(size=16,color="black"),
    legend.key=element_rect(fill=NA),
    legend.key.height=unit(1.5,"lines"),
    axis.line=element_line(color="black"),
    axis.title=element_text(size=21),
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

win.graph(80,40)
g
ggsave("Fig2.pdf")   
dev.off()
