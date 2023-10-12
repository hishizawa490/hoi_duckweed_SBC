library(tidyverse)
library(doParallel)
source("./functions.R")
data<-read_csv("data1_rawcfu.csv") %>% mutate(cfu=log10(cfu))

FLMboot<-function(i){
  set.seed(i)
  sample<-data %>% filter(nspecies<=3) %>% group_by(system,strain) %>% sample_frac(1,replace=T) %>% ungroup() 
  if(i==0){sample<-data %>% filter(nspecies <=2)}
  int<-extract_int(sample)
  preddata<-data[,1:3] %>% left_join(sample,by=c("strain","system","nspecies")) %>% 
    group_by(strain,system,nspecies) %>% summarise(cfu=mean(cfu),.groups="drop") %>% arrange(strain,nspecies,rev(system))
  for(n in 3:7){
    for(sys in unique(filter(preddata,nspecies==n)$system)){
      weight<-matrix(ncol=n,nrow=n,1)
      constant<-matrix(nrow=n,ncol=1)
      mem<-sys_to_namelist(sys)
      sub<-subsystem_all(sys,2,2)
      for(j in 1:n){
        constant[j,1]<-mean(as.matrix(filter(preddata,strain==mem[j],system %in% sub)$cfu))
        for(k in seq(1,n)[seq(1,n)!=j]){
          weight[k,j]<-(-(n-2)/(n-1))*mean(filter(int,spj==mem[j],spi==mem[k],system %in% sub)$Cij)
        }
      }
      solution<-solve(weight,constant)
      for(l in 1:n){
        preddata[preddata$strain==mem[l] & preddata$system==sys,4]<-solution[l,1]
      }
    }
  }
  preddata$cfu
}

cl<-makeCluster(12,type="PSOCK")
registerDoParallel(cl)
res<-foreach(i = 1:3000,.combine="cbind",.packages=c("tidyverse")) %dopar% (FLMboot(i))
stopCluster(cl)

sd_pred<-apply(res,MARGIN=1,sd)
pred<-FLMboot(0)

result<-data %>% group_by(strain,system,nspecies) %>% summarize(obs=mean(cfu),sd_obs=sd(cfu)) %>% 
  arrange(strain,nspecies,rev(system)) %>% ungroup() %>% mutate(pred=pred,sd_pred=sd_pred)

write.csv(result,"data7_abundance_prediction_pair.csv",row.names=FALSE)
