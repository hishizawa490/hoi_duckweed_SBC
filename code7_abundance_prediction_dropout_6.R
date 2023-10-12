library(tidyverse)
library(doParallel)
source("./functions.R")
data<-read_csv("data1_rawcfu.csv") %>% mutate(cfu=log10(cfu))
stats<-read_csv("data3_cfustats.csv")

result<-c()
for(i in 1:100){
  set.seed(i)
  for(sys in subsystem_all("s3612475",6,6)){
    sys3<-subsystem_all(sys,3,3)
  for(n_3 in 0:20){
    sample<-stats %>% filter(system %in% c(subsystem_all(sys,1,2),sample(sys3,n_3)))
    int<-extract_int(sample) %>% na.omit()
    preddata<-data[,1:3] %>% left_join(sample,by=c("strain","system","nspecies")) %>% 
      group_by(strain,system,nspecies) %>% summarise(cfu=mean(cfu),.groups="drop") %>% arrange(strain,nspecies,rev(system))
    weight<-matrix(ncol=6,nrow=6,0)
    constant<-matrix(nrow=6,ncol=1)
    mem<-sys_to_namelist(sys)
    for(j in 1:6){
      weight[j,j]<-1
      sub3<-intersect(subsystem_all(sys,3,3,mem[j]),unique(sample$system))
      sub2<-intersect(subsystem_all(sys,2,2,mem[j]),unique(sample$system))
      nodes<-c(sub3,setdiff(sub2,setdiff(subsystem_all(sub3),sub3)))
      constant[j,1]<-mean(as.matrix(filter(preddata,strain==mem[j],system %in% nodes)$cfu))
      for(k in nodes){
        subint3<-filter(int,spi==mem[j],spj %!in% sys_to_namelist(k),nspecies==3) %>% group_by(spj) %>% summarize(meanCij=mean(Cij))
        subint2<-filter(int,spi==mem[j],spj %!in% sys_to_namelist(k),nspecies==2)
        for(l in mem[-j]){
          if(l %in% subint3$spj==T){
            weight[j,which(mem==l)]<-weight[j,which(mem==l)]-filter(subint3,spj==l)$meanCij/length(nodes)
          }else if(l %in% subint2$spj ==T){
            weight[j,which(mem==l)]<-weight[j,which(mem==l)]-filter(int,spi==mem[j],spj==l,nspecies==2)$Cij/length(nodes)
          }else{
          }
        }
      }
    }
    solution<-solve(weight,constant)
    for(l in 1:6){
      preddata[preddata$strain==mem[l] & preddata$system==sys,4]<-solution[l,1]
    }
    result<-rbind(result,preddata %>% filter(system==sys)%>%mutate(n_3=n_3)%>%mutate(trial=i))
  }
  }
}
write.csv(result,"data9_abundance_prediction_sparse_6.csv",row.names=FALSE)
