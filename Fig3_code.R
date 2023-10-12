library(tidyverse)
library(RColorBrewer)
source("./functions.R")

interactions<-read_csv("data4_interactions.csv") %>% mutate(dir=str_c(spj,"→",spi))
order<-interactions %>% filter(nspecies==2) %>% arrange(Cij) %>% mutate(order=seq(1,42)) %>% select(spi,spj,order)
interactions<-left_join(interactions,order,by=c("spi","spj"))


g<-ggplot(filter(interactions,nspecies>=3))+
  geom_hline(yintercept=0,col="gray50",size=0.7,linetype=2)+
  geom_jitter(aes(x=reorder(dir,order),y=Cij,color=as.factor(nspecies)),shape=16,width=0.2,height=0,size=2.4,alpha=0.8)+
  geom_point(data=filter(interactions,nspecies==2),aes(x=dir,y=Cij),size=2.4,color="red")+
  geom_errorbar(data=filter(interactions,nspecies==2),aes(x=dir,ymin=Cij-1.96*sd_Cij,ymax=Cij+1.96*sd_Cij),width=0.2,lwd=0.8,color="red")+
  scale_color_manual(values=brewer.pal(n=6,"YlGnBu")[2:6],guide="none")+
  scale_y_continuous(name="")+
  theme_bw(base_size=20)+
  theme(
    aspect.ratio=0.27,
    axis.text.x=element_text(color="black",size=13,angle=45,vjust=1,hjust=1,margin=margin(0,0,0,0,"lines")),
    axis.text.y=element_text(color="black",size=13),
    panel.grid=element_blank(),
    legend.title=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=18)
  )
dev.off()
win.graph(80,30)
g

ggsave("Fig3.png",dpi=300)
