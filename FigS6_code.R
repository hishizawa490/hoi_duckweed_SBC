library(MLmetrics)
library(lmodel2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
setwd("C:/Users/uadgw/デスクトップ/WorkingFiles/モデル群集/syncomR/deposit")
source("./functions.R")

pred_pair<-read_csv("data7_abundance_prediction_pair.csv") %>% filter(nspecies>=4) %>% 
  mutate(min_obs=10^(obs-sd_obs)/10^6,max_obs=10^(obs+sd_obs)/10^6,min_pred=10^(pred-sd_pred)/10^6,max_pred=10^(pred+sd_pred)/10^6) %>%
  mutate(obs=10^obs/10^6,pred=10^pred/10^6)

g<-list()
for(i in 1:7){
  g[[i]]<-ggplot(data=pred_pair %>% filter(strain==stnames[i]))+
    geom_segment(aes(x=0,y=0,xend=max(c(pred,obs)),yend=max(c(pred,obs))),size=1)+
    geom_errorbar(aes(x=pred,ymin=min_obs,ymax=max_obs),size=0.7)+
    geom_errorbar(aes(xmin=min_pred,xmax=max_pred,y=obs),size=0.7)+
    geom_point(aes(x=pred,y=obs,fill=as.factor(nspecies)),shape=21,size=3.5,stroke=0.6)+
    geom_text(data=metrics(pred_pair,group="strain") %>% filter(strain==stnames[i]),hjust=-0.07,vjust=1.4,size=9,
              aes(x=-Inf,y=+Inf,label=strain))+
    geom_text(data=metrics(pred_pair,group="strain") %>% filter(strain==stnames[i]),hjust=1.1,vjust=-0.6,size=7,
              aes(x=Inf,y=-Inf,label=paste("RMSPE ",scalefun(100*rmspe),"%",sep="")))+
    scale_fill_manual(values=brewer.pal(n=6,"YlGnBu")[3:6],guide="none")+
    scale_x_continuous(name=expression(paste("Predicted abundance (",{Log[10]}," cfu ",{mg^-1},")")),labels=scalefun,breaks=breaks_fun2,expand=c(0.19,0.1))+
    scale_y_continuous(name=expression(paste("Observed abundance ( ",{Log[10]}," cfu ",{mg^-1},")")),labels=scalefun,breaks=breaks_fun2,expand=c(0.19,0.1))+
    theme_bw(base_size=30)+
    theme(
      aspect.ratio=1,
      axis.text.x=element_text(color="black",size=21,margin=margin(0.03,0.1,0.1,0.1,"lines")),
      axis.text.y=element_text(color="black",size=21,margin=margin(0.1,0.03,0.1,0.1,"lines")),
      axis.title=element_blank(),
      axis.ticks.length=unit(0.3,"lines"),
      axis.ticks=element_line(size=1),
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line=element_line(color="black",size=1.2),
      plot.margin=margin(0.2,0.2,0.2,0.2,unit="lines")
    )+
    coord_obs_pred(xlim=c(0,NA))
}

dev.off()
win.graph(20,10)
g[[1]]+g[[2]]+g[[3]]+g[[4]]+g[[5]]+g[[6]]+g[[7]]+plot_layout(ncol=4)

ggsave("FigS6.pdf")

