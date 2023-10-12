library(tidyverse)
library(cowplot)
library(reshape2)
library(agricolae)
source("./functions.R")
scalefun<-function(x){sprintf("%.1f",x)}
data<-read_csv("data1_rawcfu.csv")
for(i in 1:7){
  data<-mutate(data,!!stnames[i]:=ifelse(substr(system,i+1,i+1)!="0",1,0))
}

for(i in 1:7){
  d<-filter(data,strain==unique(data$strain)[i])
  d$cfu<-d$cfu/1000000
  d<-left_join(d,d %>% group_by(system) %>% summarise(cfuave=mean(cfu)),by="system")
  dlong<-as_tibble(melt(d,id.vars=names(d)[-6:(-12)],variable.name="vars"))

  g1<-ggplot(data=d,aes(x=reorder(system,cfuave),y=cfu))+
    geom_boxplot(lwd=0.5,fatten=0.5,outlier.shape = NA,fill=stcolors[i],color="black")+
    geom_jitter(width=0.15,height=0,size=2.3,shape=21,fill="white",stroke=0.8)+
    scale_y_continuous(name=paste(stnames[i],"colonizatoin (cfu/mg)"),labels=scalefun,limits=c(0,NA))+
    theme_bw(base_size=25)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=20,color="black"),
          panel.grid=element_blank(),
          aspect.ratio=0.32,
          plot.margin=unit(c(0.5,0.5,-0.5,0.5),"lines"))
  
  g2<-ggplot(data=filter(dlong,rep==1,vars!=stnames[i]))+
    geom_point(aes(x=reorder(system,cfuave),y="DW039"),color="transparent")+
    geom_point(aes(x=reorder(system,cfuave),y=vars,fill=paste(vars,value),color=vars),size=4.8,stroke=0.5,shape=21)+
    scale_color_manual(values=stcolors[-i],guide="none")+  
    scale_fill_manual(values=c(rep("transparent",6),stcolors[-i])[c(1,7,2,8,3,9,4,10,5,11,6,12)],guide="none")+
    scale_y_discrete(limits=rev(stnames[-i]),name="")+
    theme_void()+
    theme(axis.text.y=element_blank(),
          axis.title.y=element_text(color="black",size=16,angle=90,margin=unit(c(0,1,0.2,0),"lines")),
          aspect.ratio=0.096,
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"lines"))
  
  win.graph(75,30)
  plot_grid(g1,g2,nrow=2,align="v",rel_heights = c(0.5,0.172))
  ggsave(filename=paste("FigS3_",stnames[i],".pdf",sep=""),dpi=300)
  dev.off()
}