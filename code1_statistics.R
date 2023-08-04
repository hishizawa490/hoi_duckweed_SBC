library(vegan)
library(multcomp)
library(agricolae)
library(tidyverse)
library(fitdistrplus)
setwd("C:/Users/uadgw/デスクトップ/WorkingFiles/モデル群集/syncomR/deposit")
source("./functions.R")
data<-read_csv("data1_rawcfu.csv")

#Test for log-normal distribution
res<-tribble(~"strain",~"norm_loglik",~"norm_aic",~"lnorm_loglik",~"lnorm_aic")
for(i in stnames){
  d<-filter(data,strain==i)$cfu
  norm<-fitdist(d,"norm")
  lnorm<-fitdist(d,"lnorm")
  res<-res %>% add_row(
    "strain"=i,
    "norm_loglik"=norm$loglik,
    "norm_aic"=norm$aic,
    "lnorm_loglik"=lnorm$loglik,
    "lnorm_aic"=lnorm$aic
  )
}
write.csv(res,"data2_model_fit.csv",row.names=FALSE)


#Calculate fundamental statistics and multiple comparison test
data<-data %>% mutate(cfu=log10(cfu)) %>% mutate(system=as.factor(system))

res<-numeric()
for(i in stnames){
  d<-filter(data,strain==i) 
  tukey<-HSD.test(aov(cfu~system,data=d),"system",group=T)$groups
  tukey$system<-row.names(tukey)
  dunnett<-summary(glht(aov(cfu~system,data=d),linfct=mcp(system="Dunnett")))
  dunnett<-data.frame("system"=substr(names(dunnett$test$tstat),1,8),
                      "Dunnett_t"=dunnett$test$tstat,
                      "Dunnett_p"=dunnett$test$pvalues)
  stats<-left_join(tukey,dunnett,by="system") %>% 
    left_join(d %>% group_by(strain,nspecies,system) %>% summarize(sd=sd(cfu),n_rep=length(cfu)),by="system")
  res<-rbind(res,stats)
}
res<-res[,c(6,3,7,1,8,9,4,5,2)] %>% arrange(strain,nspecies,rev(system))
write.csv(res,"data3_cfustats.csv",row.names=FALSE)


#two-way ANOVA for the trio-based abundance prediction

pred<-read_csv("data6_abundance_prediction_trio.csv") %>% filter(nspecies>=4) %>% mutate(error=abs(10^obs-10^pred))
res<-aov(error~nspecies*strain,data=pred)
summary(res)

