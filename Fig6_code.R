library(tidyverse)
library(vegan)
source("./functions.R")

stats<-read_csv("data3_cfustats.csv") %>% mutate(obs=cfu) %>% dplyr::select(strain,system,obs) 

result<-read_csv("data8_abundance_prediction_sparse.csv") %>% gather(key="trial",value="cfu",contains("result")) %>% 
  mutate(system="s3612475",trial=sort(rep(1:100,252))) %>% left_join(stats,by=c("strain","system")) %>%
  dplyr::select(strain,system,n_3,trial,cfu,obs)

bcdis<-tribble(~"system",~"n_3",~"trial",~"bcdis")
for(sys in unique(result$system)){
  for(i in 1:100){
    for(j in 0:35){
      res<-result %>% filter(system==sys,trial==i,n_3==j)
      bcdis<-bcdis %>% add_row(
        "system"=sys,
        "n_3"=j,
        "trial"=i,
        "bcdis"=as.numeric(vegdist(rbind(10^res$cfu,10^res$obs)))
      )
    }
  }
}

g1<-ggplot(data=na.omit(bcdis),aes(x=n_3,y=bcdis))+
  geom_jitter(alpha=0.15,size=2.7,height=0,width=0.18,col="cyan4")+
  stat_summary(fun=median,geom="line",aes(group=1),col="red",size=1.8)+
  scale_x_continuous(name="N of trio combinations",expand=c(0.02,0.02))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6),name="Bray-Curtis dissimilarity",expand=c(0.004,0),limits=c(0,0.75))+
  theme_bw(base_size=26)+
  theme(
    aspect.ratio=0.5,
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line=element_line(color="black",size=1.2),
    axis.title=element_text(size=18),
    axis.text=element_text(color="black",size=18),
    strip.text=element_text(size=17,margin=margin(0,0,0,0,"lines")),
    strip.background=element_blank(),
    plot.margin=unit(c(0.5,0.5,-0.2,0.5),"lines")
  )

sys3<-subsystem_all("s3612475",3,3)
data<-read_csv("data10_prediction_improvement.csv") 

g2<-ggplot(data=data)+
  geom_boxplot(aes(x=reorder(sys,bcdis-wo_bcdis),y=wo_bcdis-bcdis),lwd=0.4,fatten=1,fill="cyan4",alpha=0.8,outlier.size=0.6)+
  scale_y_continuous(name="Prediction Improvement")+
  scale_x_discrete(expand=c(0.02,0))+
  theme_bw(base_size=26)+
  theme(
    panel.border=element_blank(),
    axis.line=element_line(color="black",size=1.2),axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title=element_text(size=18),
    axis.text.y=element_text(size=18,color="black"),
    panel.grid=element_blank(),
    aspect.ratio=0.5,
    plot.margin=unit(c(0,0.5,-0.5,0.5),"lines")
    )

pdata<-data %>% group_by(sys) %>% summarize(imp=mean(bcdis-wo_bcdis)) %>%
  rowwise() %>% mutate(vars=list(sys_to_namelist(sys))) %>% unnest(vars)

g3<-ggplot(data=pdata)+
  geom_point(aes(x=reorder(sys,imp),y="DW039"),color="transparent")+
  geom_point(aes(x=reorder(sys,imp),y=vars,fill=vars,color=vars),size=2.5,stroke=0.5,shape=21)+
  scale_x_discrete(name="Composition of trio combinations",expand=c(0.02,0))+
  scale_color_manual(values=stcolors,guide="none")+  
  scale_fill_manual(values=stcolors)+
  scale_y_discrete(limits=rev(stnames),name="")+
  guides(fill=guide_legend(ncol=2))+
  theme_void()+
  theme(axis.text.y=element_blank(),
        legend.position=c(0.8,4.2),
        legend.title=element_blank(),
        axis.title.x=element_text(size=18,color="black",vjust=-0.1),
        legend.text=element_text(size=13,color="black"),
        aspect.ratio=0.13,
        plot.margin=unit(c(-0.7,0.5,0.5,0.5),"lines"))

dev.off()
win.graph(10,12)
plot_grid(g1,g2,g3,nrow=3,align="v",rel_heights = c(0.6,0.5,0.18))

ggsave("Fig6.pdf")