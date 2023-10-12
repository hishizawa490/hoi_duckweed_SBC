library(tidyverse)
library(vegan)
source("./functions.R")

data_trio<-read_csv("data6_abundance_prediction_trio.csv")
data_pair<-read_csv("data7_abundance_prediction_pair.csv")
data<-read_csv("data3_cfustats.csv") 

res<-tribble(~"system",~"nspecies",~"bcdis_trio",~"bcdis_pair")
for(i in unique(filter(data,nspecies>=4)$system)){
  pred_trio<-10^filter(data_trio,system==i)$pred 
  pred_pair<-10^filter(data_pair,system==i)$pred
  obs<-10^filter(data,system==i)$cfu
 
  res<-res %>% add_row(
    "system"=i,
    "nspecies"=sys_to_n(i),
    "bcdis_trio"=as.numeric(vegdist(rbind(pred_trio,obs))),
    "bcdis_pair"=as.numeric(vegdist(rbind(pred_pair,obs)))
  )
}

res<-pivot_longer(res,cols=c("bcdis_trio","bcdis_pair"),values_to="bcdis",names_to="model")
group<-HSD.test(aov(bcdis~sys,data=mutate(res,sys=str_c(nspecies,model))),"sys",group=T)$groups
tukey<-res %>% group_by(nspecies,model) %>% reframe(max=max(bcdis)) %>% 
  mutate(sys=str_c(nspecies,model)) %>% left_join(rownames_to_column(group,"sys"),"sys")

g<-ggplot(data=res)+
  geom_boxplot(aes(x=as.factor(nspecies),y=bcdis,fill=model),position=position_dodge(1),alpha=0.8,size=1)+
  geom_jitter(aes(x=as.factor(nspecies),y=bcdis,fill=model),position=position_dodge(1),shape=21,size=3,stroke=1)+
  geom_text(data=filter(tukey,model=="bcdis_pair"),aes(x=as.factor(nspecies),y=max+0.05,label=groups),size=9,nudge_x=-0.25)+
  geom_text(data=filter(tukey,model=="bcdis_trio"),aes(x=as.factor(nspecies),y=max+0.05,label=groups),size=9,nudge_x=0.25)+
  scale_fill_manual(values=c("hotpink2","lightseagreen"),guide=NULL)+
  scale_y_continuous(limits=c(0,0.79),expand=c(0,0))+
  theme_bw(base_size=34)+
  theme(
    aspect.ratio=1,
    panel.grid=element_blank(),
    axis.text=element_text(color="black",size=32),
    axis.title=element_blank(),
    panel.border=element_blank(),
    axis.line=element_line(color="black",size=2)
  )
dev.off()
win.graph(10,10)
g
ggsave("Fig.5d.pdf")