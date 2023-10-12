library(tidyverse)
library(reshape2)
source("./functions.R")
scalefun<-function(x){sprintf("%.1f",x)}
data<-read_csv("data10_biased_inoculation.csv")

data$cfu<-data$cfu/1000000
d<-left_join(data,data %>% group_by(strain,system) %>% summarise(cfuave=mean(cfu),sd=sd(cfu)),by=c("strain","system"))  %>% 
  filter(rep==1) %>% arrange(desc(system),strain) %>% dplyr::select(-rep,-cfu) %>% mutate("order"=1,"cfusum"=cfuave)
d$strain<-factor(d$strain,levels=rev(stnames))

order<-1; val<-d$cfuave[1]
for(i in 2:nrow(d)){
  if(d$system[i]!=d$system[i-1]){
    val<-d$cfuave[i]
    order<-order+1
  }else{val<-val+d$cfuave[i]}
  d$cfusum[i]<-val
  d$order[i]<-order
}

g<-ggplot(data=d)+
  geom_bar(aes(x=system,y=cfuave),fill="transparent",stat="identity")+
  geom_bar(aes(x=system,y=cfuave,fill=strain),color="black",stat="identity",lwd=0.01)+
  geom_errorbar(aes(x=system,ymin=cfusum-sd,ymax=cfusum+sd),color="grey10",width=0.2,lwd=0.3)+
  scale_y_continuous(labels=scalefun,limits=c(NA,NA),expand=c(0,0))+
  scale_fill_manual(values=rev(stcolors))+
  theme_bw(base_size=25)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(color="black",size=19),
        panel.grid=element_blank(),
        legend.title=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(color="black",size=1),
        aspect.ratio=1,
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"lines"))+
  coord_cartesian(ylim=c(0,6.2),xlim=c(1,7))
  
win.graph(10,10)
g
  
ggsave(filename="FigS2.pdf")
dev.off()
