library(tidyverse)
library(agricolae)
library(PTXQC)
source("./functions.R")
data<-read_csv("data1_rawcfu.csv") %>% mutate(cfu=log10(cfu))

interactions<-tribble(~"spi",~"spj",~"system",~"nspecies",~"Pi.j",~"Pi",~"Pj",~"Cij",~"sd_Cij",~"hsd")
for(i in stnames){
  d<-filter(data,strain==i)
  tukey<-HSD.test(aov(cfu~system,data=d),"system",group=T)$groups
  tukey$system<-row.names(tukey)
  d<-left_join(d,tukey[,c(3,2)],by="system")
  for(sys in unique(filter(d,nspecies>=2)$system)){
    sublist<-subsystem(sys)
    for(sub in filter(sublist,add!=i)$subsystem){
      Pi.j<-filter(d,system==sys)$cfu
      Pi<-filter(d,system==sub)$cfu
      Pj<-filter(data,system==sys,strain==filter(sublist,subsystem==sub)$add)$cfu
      coefs<-numeric(3000)
      for(k in 1:3000){
        set.seed(k)
        coefs[k]<-(mean(sample(Pi.j,length(Pi.j),replace=T))-mean(sample(Pi,length(Pi),replace=T)))/mean(sample(Pj,length(Pj),replace=T))
      }
      interactions<-interactions %>% 
        add_row("spi"=i,
                "spj"=filter(sublist,subsystem==sub)$add,
                "system"=sys,
                "nspecies"=sys_to_n(sys),
                "Pi.j"=mean(Pi.j),
                "Pi"=mean(Pi),
                "Pj"=mean(Pj),
                "Cij"=(mean(Pi.j)-mean(Pi))/mean(Pj),
                "sd_Cij"=sd(coefs),
                "hsd"=ifelse(LCSn(c(filter(d,system==sys)$groups[1],filter(d,system==sub)$groups[1]))=="","sig","n.s."))
    }
  }
}

write.csv(interactions,"data4_interactions.csv",row.names=FALSE)

interactions<-interactions %>% left_join(interactions%>%filter(nspecies==2)%>%select(spi,spj,Cij),by=c("spi","spj")) %>% 
  mutate(dCij=Cij.x-Cij.y) %>% filter(nspecies >=3)
res<-aov(dCij~dir,data=interactions)
summary(res)
