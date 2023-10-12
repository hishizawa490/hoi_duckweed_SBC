stnames<-c("DW039","DW067","DW100","DW102","DW145","DW147","DW155")
stcolors<-c("#69C232","#FFC000","#9148C8","#FA8738","#F84447","#35978E","#323FB8")
strainlist<-c("DW100","DW102","DW039","DW145","DW155","DW067","DW147")
fundniche<-c(3600989.1,515835.4,311853.4,1409510.3,1694884.2,713002.6,1378926.5)
fundniche_sd<-c(445200.85,78789.57,21362.76,297663.98,130778.35,114255.09,185237.51)
'%!in%'=Negate('%in%')

scalefun<-function(x){sprintf("%.1f",x)}

breaks_fun<-function(y){
  if((max(y)-min(y))<1.2){
    seq(2*round(min(y/2),1),round(max(y),1),by=0.2)
  }else{
    seq(5*round(min(y/5),1),round(max(y),1),by=0.5)
  }
}

breaks_fun2<-function(y){
  seq<-seq(0,round(max(y),1)+0.2,by=0.2)
  if(length(seq)>6){seq<-seq(0,round(max(y),1)+0.5,by=0.5)}
  if(length(seq)>6){seq<-seq(0,round(max(y),1)+1,by=1)}
  if(length(seq)>6){seq<-seq(0,round(max(y),1)+2,by=2)}
  if(length(seq)>6){seq<-seq(0,round(max(y),1)+3,by=3)}
  seq
}

#Output subsystem ID and omitted species name
subsystem<-function(system,n=1){
  subs<-numeric();subs2<-numeric()
  for(i in c(3,6,1,2,4,7,5)){
   sub<-gsub(i,"0",system) 
   if(sub!=system && sub!="s0000000"){
     subs<-rbind(subs,c(sub,strainlist[i]))}
  }
  subs<-data.frame(subs)
  colnames(subs)<-c("subsystem","add")
  if(n==2){
    for(j in 1:nrow(subs)){
      for(i in c(3,6,1,2,4,7,5)[name_to_num(subs$add[j]):7]){
        sub<-gsub(i,"0",subs$subsystem[j]) 
        if(sub!=subs$subsystem[j] && sub!="s0000000"){
          subs2<-rbind(subs2,c(sub,subs$add[j],strainlist[i]))
        }
      }
      colnames(subs2)<-c("subsystem","add1","add2")
    }
  data.frame(subs2)
  }else{data.frame(subs)}
}

#Output all subsystem IDs
subsystem_all<-function(system,min=1,max=7,species=""){
  subs<-system
  for(i in c(3,6,1,2,4,7,5)){
    subs<-c(subs,gsub(i,0,subs))
  }
  subs<-unique(subs)
  subs<-subs[subs!="s0000000"]
  subs<-subs[between(sys_to_n(subs),min,max)==T]
  if(species !=""){
    for(i in 1:length(subs)){
      if(substr(subs[i],name_to_num(species)+1,name_to_num(species)+1)==0){
        subs[i]<-NA
      }
    }
  }
  na.omit(subs)
}

#calculate prediction performance
metrics<-function(data,min=4,max=7,group=""){
  data<-data %>% filter(nspecies>=min,nspecies<=max)  
  if(group[1]!=""){data<-data%>%group_by_at(group)}
  data %>% 
  summarize(rmse=RMSE(pred,obs),
            rmspe=RMSPE(pred,obs),
            mae=MAE(pred,obs),
            mape=100*MAPE(pred,obs),
            pearsonr=cor(pred,obs),
            ma_slope=tryCatch({lmodel2(obs~pred)$regression.results[2,3]},error=function(e){NA}),
            ols_slope=tryCatch({lmodel2(obs~pred)$regression.results[1,3]},error=function(e){NA}),
            n=length(pred),
            meanobs=mean(abs(obs)),
            meanpred=mean(abs(pred)))
}


#Convert cfu count table into interaction table
extract_int<-function(data){
  intlist<-tribble(~"spi",~"spj",~"system",~"subsystem",~"nspecies",~"Pi.j",~"Pi",~"Pj",~"Cij")
  col<-na.omit(match(c("cfu","pred"),colnames(data)))[1]
  for(i in unique(data$strain)){
    d<-filter(data,strain==i)
    for(sys in unique(filter(d,nspecies>=2)$system)){
      sublist<-subsystem(sys)
      for(sub in filter(sublist,add!=i)$subsystem){
        Pi.j<-filter(d,system==sys)[,col];if(is.na(sum(Pi.j))==T){Pi.j<-filter(d,system==sys)[,col-2]}
        Pi<-filter(d,system==sub)[,col];if(is.na(sum(Pi))==T){Pi<-filter(d,system==sub)[,col-2]}
        Pj<-filter(data,system==sys,strain==filter(sublist,subsystem==sub)$add)[,col]
          if(is.na(sum(Pj))==T){Pj<-filter(data,system==sys,strain==filter(sublist,subsystem==sub)$add)[,col-2]}
          intlist<-intlist %>% 
          add_row("spi"=i,
                  "spj"=filter(sublist,subsystem==sub)[1,2],
                  "system"=sys,
                  "subsystem"=sub,
                  "nspecies"=sys_to_n(sys),
                  "Pi.j"=mean(as.matrix(Pi.j)),
                  "Pi"=mean(as.matrix(Pi)),
                  "Pj"=mean(as.matrix(Pj)),
                  "Cij"=(mean(as.matrix(Pi.j))-mean(as.matrix(Pi)))/mean(as.matrix(Pj)))
      }
    }
  }
intlist %>% arrange(spi,nspecies,rev(system))
}

#Convert system ID into number of species
sys_to_n<-function(system){
  n<-7-str_count(system,pattern="0")
  n
}

#Convert system ID into the list of existing species
sys_to_namelist<-function(system){
  x<-numeric()
  for(i in 2:8){
    num<-as.numeric(substr(system,i,i))
    if(num>0){
      x<-c(x,strainlist[num])
    }
  }
  x
}

#Convert system ID into the list of co-existing species
sys_to_coex<-function(system,strain){
  x<-numeric()
  for(i in 2:8){
    num<-as.numeric(substr(system,i,i))
    if(num>0){
      if(strainlist[num]!=strain){
        x<-c(x,strainlist[num])
      }
    }
  }
  x
}

#Convert species name into numerical ID
name_to_num<-function(name){
  x<-numeric(length(name))
  for(i in 1:length(name)){
    x[i]<-which(stnames==name[i])
  }
  x
}

#Convert species name list into system ID
name_to_sys<-function(names){
  sys<-"s"
  for(i in 1:7){
    if(length(grep(stnames[i],names))>0){
      sys<-str_c(sys,grep(stnames[i],strainlist))
    }else{
      sys<-str_c(sys,"0")
    }
  }
  sys
}





scientific_notation<-function(x) {
  x <- sprintf('%.1e',x)
  x <- gsub("^(.*)e", "'\\1'e",x)
  x <- gsub("e", "%*%10^", x)
  x <- gsub('\\+', '', x)
  parse(text = x)}



#
#https://github.com/tidymodels/tune/blob/main/R/coord_obs_pred.R
CoordObsPred <-
  ggplot2::ggproto(
    "CoordObsPred",
    ggplot2::CoordFixed,
    setup_panel_params = function(self, scale_x, scale_y, params = list()) {
      # coord limits take precedence over scale limits
      rngs <- range(self$limits$x %||% scale_x$get_limits(),
                    self$limits$y %||% scale_y$get_limits(),
                    na.rm = TRUE)
      self$limits$y <- rngs
      self$limits$x <- rngs
      ggplot2::ggproto_parent(CoordFixed, self)$setup_panel_params(scale_x, scale_y, params)
    },
    aspect = function(self, ranges) {
      1 / self$ratio
    }
  )


coord_obs_pred <-
  function(ratio = 1,
           xlim = NULL,
           ylim = NULL,
           expand = TRUE,
           clip = "on") {
    ggplot2::ggproto(
      NULL,
      CoordObsPred,
      limits = list(x = xlim, y = ylim),
      ratio = ratio,
      expand = expand,
      clip = clip
    )
  }


