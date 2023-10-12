library(ggpmisc)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(lmodel2)
library(ggpmisc)
cpred<-read_csv("data5_coef_prediction.csv")

g<-list()
k<-0
for(i in 2:6){
  for(j in seq(i+1,7)){
    k<-k+1
    pdata<-cpred %>% filter(nspecies==j,Cij_from==i)
    g[[k]]<-ggplot(data=pdata,aes(x=pred,y=obs))+
      geom_point(size=1.2,col=brewer.pal(n=6,"YlGnBu")[j-1],alpha=0.8)+
      stat_ma_line(col="deeppink",size=1.5)+
      annotate(geom="text",x=Inf,y=-Inf,label=format(round(lmodel2(formula=pdata$obs~pdata$pred)$regression.results$Slope[2],2),nsmall=2),
               vjust=-1.3,hjust=1.2,size=6,col="deeppink")+
      theme_void(base_size=8)+
      theme(aspect.ratio=1,
            strip.background=element_blank(),
            axis.title=element_blank(),
            panel.grid=element_blank(),
            axis.ticks=element_blank(),
            plot.background=element_rect(fill="gray95",color="white",size=1),
            plot.margin=margin(0,0,0,0,"lines"),
            axis.text=element_blank())+
      coord_cartesian(xlim=c(-0.19,0.12),ylim=c(-0.19,0.12))
  }
}

g[[16]]<-ggplot(data=tibble("x"=1:5,"y"=1:5),aes(x=x,y=y))+
          theme_void()+
          scale_x_continuous(expand=c(0,0))+
          scale_y_continuous(expand=c(0,0))+
          theme(
            panel.background=element_rect(fill="gray80",color="white",size=2),
            panel.grid=element_blank(),
            axis.ticks=element_line(color="black",size=1),
            axis.ticks.length=unit(0.2,"lines"),
            plot.margin=margin(1,1,1,1,"lines")
          )+
          coord_cartesian(xlim=c(0.5,5.5),ylim=c(0.5,5.5))

dev.off()
win.graph(20,20)

g[[16]]+inset_element(g[[5]],0.0,0.8,0.2,1.0)+inset_element(g[[9]],0.2,0.8,0.4,1.0)+inset_element(g[[12]],0.4,0.8,0.6,1.0)+inset_element(g[[14]],0.6,0.8,0.8,1.0)+inset_element(g[[15]],0.8,0.8,1.0,1.0)+
        inset_element(g[[4]],0.0,0.6,0.2,0.8)+inset_element(g[[8]],0.2,0.6,0.4,0.8)+inset_element(g[[11]],0.4,0.6,0.6,0.8)+inset_element(g[[13]],0.6,0.6,0.8,0.8)+
        inset_element(g[[3]],0.0,0.4,0.2,0.6)+inset_element(g[[7]],0.2,0.4,0.4,0.6)+inset_element(g[[10]],0.4,0.4,0.6,0.6)+
        inset_element(g[[2]],0.0,0.2,0.2,0.4)+inset_element(g[[6]],0.2,0.2,0.4,0.4)+
        inset_element(g[[1]],0.0,0.0,0.2,0.2)

ggsave("Fig4b.pdf")