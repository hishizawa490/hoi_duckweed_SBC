library(tidyverse)
source("./functions.R")
int<-read_csv("data4_interactions.csv")
pred3<-read_csv("data6_abundance_prediction_trio.csv")
predint3<-extract_int(pred3)

int[,5:7]<-10^int[5:7]/1000000
predint3[,6:8]<-10^predint3[6:8]/1000000
max<-left_join(int %>% group_by(spi) %>% summarize(maxobs=max(Pi.j)),predint3 %>% group_by(spi) %>% summarize(maxpred=max(Pi.j))) %>% 
  rowwise() %>% mutate(max=max(maxobs,maxpred)) %>% filter(spi%in%c("DW067","DW145"))

predint3<-filter(predint3,spi%in%c("DW067","DW145"))

g<-ggplot(data=predint3)+
  geom_segment(aes(x=nspecies-1,xend=nspecies,y=Pi,yend=Pi.j,color=str_c("+ ",spj),linetype=as.factor(nspecies),size=as.factor(nspecies)))+
  geom_point(data=filter(predint3,nspecies>=4),aes(x=nspecies,y=Pi.j),size=2,shape=21,fill="white",stroke=1,color="gray40")+
  geom_point(data=filter(predint3,nspecies<=3),aes(x=nspecies,y=Pi.j),size=2,shape=21,fill="white",stroke=1.2)+
  geom_point(data=filter(predint3,nspecies==2),aes(x=1,y=Pi),size=2,shape=21,fill="white",stroke=0.8)+
  geom_point(data=filter(max,spi==spi),aes(x=1,y=max*1.08),fill="transparent",color="transparent")+
  scale_color_manual(values=stcolors,name="",guide="none")+
  scale_x_continuous(breaks=seq(1,max(int$nspecies),1))+
  scale_y_continuous(labels=scalefun,limits=c(0,NA),breaks=breaks_fun2,expand=c(0,0))+
  scale_linetype_manual(values=c(1,1,5,5,5,5,5),guide="none")+
  scale_size_manual(values=c(0.9,0.9,0.5,0.5,0.5,0.5,0.5),guide="none")+
  theme_bw(base_size=22)+
  theme(
    axis.title=element_blank(),
    axis.text.x=element_text(size=14,color="black",margin=margin(0.1,0.1,0.1,0.1,"lines")),
    axis.text.y=element_text(size=14,color="black",margin=margin(0.1,0.1,0.1,0.1,"lines")),
    axis.ticks.length=unit(0.2,"lines"),
    axis.line=element_line(color="black",size=0.7),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    aspect.ratio=0.85,
    strip.text=element_text(size=17,margin=margin(0,0,0,0,"lines")),
    strip.background=element_blank(),
    panel.spacing.x=unit(0.3,"lines"),
    panel.spacing.y=unit(0.3,"lines")
  )+
  coord_cartesian(ylim=c(-0.015,NA))+
  facet_wrap(~spi,nrow=1,scales="free")

  win.graph(10,10)
  g
  ggsave("Fig5b.pdf")
  dev.off()

