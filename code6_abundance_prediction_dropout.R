library(tidyverse)
library(doParallel)
setwd("C:/Users/uadgw/デスクトップ/WorkingFiles/モデル群集/syncomR/deposit")
source("./functions.R")
data<-read_csv("data1_rawcfu.csv") %>% mutate(cfu=log10(cfu))
stats<-read_csv("data3_cfustats.csv")

sys3<-subsystem_all("s3612475",3,3)

boot<-function(i){  
  result<-c()
  for(n_3 in 0:35){
    set.seed(i)
    sample<-stats %>% filter(system %in% c(subsystem_all("s3612475",1,2),sample(sys3,n_3)))
    int<-extract_int(sample) %>% na.omit()
    preddata<-data[,1:3] %>% left_join(sample,by=c("strain","system","nspecies")) %>% 
      group_by(strain,system,nspecies) %>% summarise(cfu=mean(cfu),.groups="drop") %>% arrange(strain,nspecies,rev(system))
    weight<-matrix(ncol=7,nrow=7,0)
    constant<-matrix(nrow=7,ncol=1)
    mem<-sys_to_namelist("S3612475")
    for(j in 1:7){
      weight[j,j]<-1
      sub3<-intersect(subsystem_all("s3612475",3,3,mem[j]),unique(sample$system))
      sub2<-intersect(subsystem_all("s3612475",2,2,mem[j]),unique(sample$system))
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
    for(l in 1:7){
      preddata[preddata$strain==mem[l] & preddata$system=="s3612475",4]<-solution[l,1]
    }
    result<-rbind(result,preddata %>% filter(nspecies==7))
  }
  result$cfu
}


cl<-makeCluster(12,type="PSOCK")
registerDoParallel(cl)
res<-foreach(i = 1:100,.combine="cbind",.packages=c("tidyverse")) %dopar% (boot(i))
stopCluster(cl)

sd_pred<-apply(res,MARGIN=1,sd)

result<-tibble("strain"=rep(stnames,36),"n_3"=sort(rep(0:35,7)))
result<-as_tibble(cbind(result,res,sd_pred))

write.csv(result,"data8_abundance_prediction_sparse.csv",row.names=FALSE)


