library(tidyverse)
library(reshape2)
source("./functions.R")

pred_tri<-read_csv("data6_abundance_prediction_trio.csv") %>% filter(nspecies>=4)
pred_bin<-read_csv("data7_abundance_prediction_pair.csv")

bars<-tibble("Observed"=10^filter(pred_tri,nspecies==7)$obs/1000000,
             "Pred_trio"=10^filter(pred_tri,nspecies==7)$pred/1000000,
             "Pred_pair"=10^filter(pred_bin,nspecies==7)$pred/1000000,
             "Null"=fundniche/1000000,
             "strain"=stnames) %>%
  melt() %>% mutate(variable=factor(variable,levels=c("Null","Pred_pair","Pred_trio","Observed")),strain=factor(strain,levels=rev(stnames)))

g<-ggplot(data=bars)+
  geom_bar(aes(x=value,y=variable,fill=strain),stat="identity",color="black",lwd=1,width=0.8)+
  scale_x_continuous(expand=c(0,0),limits=c(0,14.2))+
  scale_fill_manual(values=rev(stcolors),name="")+
  guides(fill=guide_legend(ncol=3,reverse=TRUE))+
  theme_bw(base_size=40)+
  theme(
    aspect.ratio=0.56,
    axis.title=element_blank(),
    axis.text.x=element_text(color="black",size=24,margin=margin(0.1,0.1,0.1,0.1,"lines")),
    axis.text.y=element_text(color="black",size=24,margin=margin(0.1,0.1,0.1,0.1,"lines")),
    legend.position=c(0.72,0.9),
    legend.text=element_text(size=20,margin=margin(0,0,0.1,-1,"lines")),
    legend.title=element_blank(),
    legend.background=element_rect(fill="transparent"),
    panel.grid=element_blank(),
    axis.ticks.length=unit(0.3,"lines"),
    panel.border=element_blank(),
    axis.line=element_line(color="black",size=1.5)
  )
win.graph(10,7)
g

ggsave("Fig5e.pdf")