######################
### Function drawing graph with stem excentricity in each level
###	
### 
######################

drawExcentricityGraph<-function(tree,plot,heights,trw,extc,method="Schweingruber"){
  Schw<-exc$Schweingruber
  Braa<-exc$Braam
  Ales<-exc$Alestalo
  
  if(method=="Schweingruber"){
    g1<-.plotGraph("1","1",heights,trw,Schw,direction="North-South")
    g2<-.plotGraph("1","1",heights,trw,Schw,direction="East-West")
    lab<-c("Schweingruber - N-S","Schweingruber W-E")
    #plot_grid(g1, g2, labels=lab, ncol = 2, nrow = 1)
  }else{;}
  if(method=="Braam"){
    g1<-.plotGraph("1","1",heights,trw,Braa,direction="North-South")
    g2<-.plotGraph("1","1",heights,trw,Braa,direction="East-West")
    lab<-c("Braam - N-S","Braam W-E")
  }else{;}
  if(method=="Alestalo"){
    g1<-.plotGraph("1","1",heights,trw,Ales,direction="North-South")
    g2<-.plotGraph("1","1",heights,trw,Ales,direction="East-West")
    lab<-c("Alestalo - N-S","Alestalo W-E")
  }else{;}
  plot_grid(g1, g2, labels=lab, ncol = 2, nrow = 1)
}

######################
### Function drawing graph of cross-section profile
###	
### 
######################

drawCrossSectionProfile<-function(trw,tree=1,plot=1,level=1){
  w <- .seriesTRW_two(trw,tree,plot,level)
  w <- .replace(w)
  widths<-.widthsCalculation(w)
  widths$S<-widths$S*-1;widths$E<-widths$E*-1
  elipse<-.elipseFun(widths[1,])
  for(i in 2:length(widths$cambAge)){
    elipse2<-.elipseFun(widths[i,])
    elipse<-rbind(elipse,elipse2)
  }
  graph<-ggplot(elipse,aes(x,y)) + geom_point(aes(colour = factor(calYear)))+geom_path(aes(colour = factor(calYear)),size=1.2)
  graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  graph
  return(graph)
}

######################
### Function drawing graph of stem alometry
###	
### 
######################

drawAlometryGraph<-function(trw,vyska,plot=1,tree=1,dir="N-S"){
  d<-.dataGraphAlometry(trw,vysky,plot,tree,dir)
  d <- d[order(d$curve,d$e),] 
  graph<-ggplot(d,aes(trwSum,abs(height))) + geom_point(aes(colour = factor(year)),size=3)+geom_path(aes(colour = factor(year)),size=1.2)
  graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),  axis.line = element_line(colour = "black", size = 1, linetype = "solid"))
  graph
  return(graph)
}

