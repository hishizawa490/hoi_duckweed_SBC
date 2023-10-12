library(tidyverse)
library(lmodel2)
library(MLmetrics)
source("./functions.R")
setwd("C:/Users/uadgw/デスクトップ/WorkingFiles/モデル群集/syncomR/deposit")

cpred<-read_csv("data5_coef_prediction.csv")
score<-metrics(cpred,min=3,max=7,group=c("Cij_from","nspecies")) %>% ungroup() 

g<-ggplot(data=score)+
  geom_tile(aes(x=as.factor(Cij_from),y=as.factor(nspecies),fill=rmse),color="white",lwd=1)+
  geom_text(aes(x=as.factor(Cij_from),y=as.factor(nspecies),label=format(round(rmse,3),nsmall=3)),size=8,vjust=0.5)+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0),breaks=3:7)+
  scale_fill_gradient(low="white",high="#FF4444",limits=c(0.015,0.085),guide=NULL)+
  theme_minimal(base_size=20)+
  theme(
    aspect.ratio=1,
    axis.title=element_blank(),
    axis.text=element_blank(),
    panel.grid=element_blank(),
    panel.background=element_rect(fill="gray80",color="white"),
    axis.ticks=element_line(color="black")
  )
dev.off()
win.graph(20,20)
g

ggsave("Fig4d_coef_pred_rmse.pdf")
