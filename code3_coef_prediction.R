library(tidyverse)
source("./functions.R")

int<-read_csv("data4_interactions.csv")

cpred<-tribble(~"spi",~"spj",~"system",~"nspecies",~"Cij_from",~"obs",~"sd_obs",~"pred",~"sd_pred",~"n_Cij_from")
for(n in 2:6){
  for(i in 1:nrow(filter(int,nspecies>=n+1))){
      topred<-filter(int,nspecies>=n+1)[i,]
      subint<-filter(int,spi==topred$spi,spj==topred$spj,system %in% subsystem_all(topred$system,n,n))
      cpred<-cpred %>% add_row(
        "spi"=topred$spi,
        "spj"=topred$spj,
        "system"=topred$system,
        "Cij_from"=n,
        "nspecies"=topred$nspecies,
        "obs"=topred$Cij,
        "sd_obs"=topred$sd_Cij,
        "pred"=mean(subint$Cij),
        "sd_pred"=sqrt(sum(subint$sd_Cij^2)),
        "n_Cij_from"=length(subint$Cij)
      )
  }
}

write.csv(cpred,"data5_coef_prediction.csv",row.names=F)
