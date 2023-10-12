library(tidyverse)
library(vegan)
source("./functions.R")

stats<-read_csv("data3_cfustats.csv") %>% mutate(obs=cfu) %>% dplyr::select(strain,system,obs) 

result<-read_csv("data9_abundance_prediction_sparse_6.csv") %>% left_join(stats,by=c("strain","system")) %>%
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

pdata<-bcdis %>% rowwise() %>% mutate(system=str_c("-",setdiff(stnames,(sys_to_namelist(system)))))

g<-ggplot(data=na.omit(pdata),aes(x=n_3,y=bcdis))+
  geom_jitter(alpha=0.15,size=1.7,height=0,width=0.18,col="cyan4")+
  stat_summary(fun=median,geom="line",aes(group=1),col="red",size=1.4)+
  scale_x_continuous(name="N of trio combinations",expand=c(0.02,0.02))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6),name="Bray-Curtis dissimilarity",expand=c(0.004,0),limits=c(0,0.75))+
  theme_bw(base_size=22)+
  theme(
    aspect.ratio=0.7,
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line=element_line(color="black",size=1.2),
    axis.text=element_text(color="black",size=14),
    strip.text=element_text(size=17,margin=margin(0,0,0,0,"lines")),
    strip.background=element_blank()
  )+
  facet_wrap(nrow=2,scales="free",~factor(system,levels=c("-DW039","-DW067","-DW100","-DW102","-DW145","-DW147","-DW155")))
dev.off()
win.graph(20,10)
g
ggsave("FigS7.pdf")