library(tidyverse)
library(cowplot)
library(reshape2)
library(agricolae)
source("./functions.R")
scalefun<-function(x){sprintf("%.1f",x)}
data<-read_csv("data1_rawcfu.csv")

data$cfu<-data$cfu/1000000
d<-left_join(data,data %>% group_by(strain,system) %>% summarise(cfuave=mean(cfu),sd=sd(cfu)),by=c("strain","system"))  %>% 
  filter(rep==1) %>% arrange(nspecies,desc(system),strain) %>% dplyr::select(-rep,-cfu) %>% mutate("order"=1,"cfusum"=cfuave)
d$strain<-factor(d$strain,levels=rev(stnames))

order<-1
for(i in 2:nrow(d)){
  if(d$system[i]!=d$system[i-1]){
    val<-d$cfuave[i]
    order<-order+1
  }else{val<-val+d$cfuave[i]}
  d$cfusum[i]<-val
  d$order[i]<-order
}
  g<-ggplot(data=d)+
    geom_bar(aes(x=order,y=cfuave),fill="transparent",stat="identity")+
    geom_rect(aes(xmin=0.5,xmax=7.5,ymin=0,ymax=7*10^6),fill="gray93",alpha=0.1)+
    geom_rect(aes(xmin=28.5,xmax=63.5,ymin=0,ymax=7*10^6),fill="gray93",alpha=0.1)+
    geom_rect(aes(xmin=98.5,xmax=119.5,ymin=0,ymax=7*10^6),fill="gray93",alpha=0.1)+
    geom_rect(aes(xmin=126.5,xmax=127.5,ymin=0,ymax=7*10^6),fill="gray93",alpha=0.1)+
    geom_bar(aes(x=order,y=cfuave,fill=strain),color="white",stat="identity",lwd=0.01)+
    geom_errorbar(aes(x=order,ymin=cfusum-sd,ymax=cfusum+sd),color="grey10",width=0.5,lwd=0.3)+
    scale_y_continuous(labels=scalefun,limits=c(NA,NA))+
    scale_x_continuous(breaks=c(0.5,7.5,28.5,63.5,98.5,119.5,126.5,127.5))+
    scale_fill_manual(values=rev(stcolors),guide="none")+
    theme_bw(base_size=25)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(color="black",size=19),
          panel.grid=element_blank(),
          aspect.ratio=0.32,
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"lines"))+
    coord_cartesian(ylim=c(0.3,6.7),xlim=c(6,122))
  
  win.graph(75,23.5)
  g
  
  ggsave("FigS3_total.pdf")
  dev.off()
