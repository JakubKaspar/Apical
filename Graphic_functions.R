######################
### Function drawing graph with stem excentricity in each level
###	
### 
######################

drawExcentricityGraph<-function(plot,tree,heights,trw,exc,withAlometry=T,method="Schweingruber"){
  Schw<-exc$Schweingruber
  Braa<-exc$Braam
  Ales<-exc$Alestalo
  
  heights<-subset(heights,Plocha_ID==plot & Strom_ID==tree)
  height<-heights$Vyska_cm
  #print(heights)
  
  m<-subset(.IDdistinct(trw)$IDAspect,IDPlot==plot & IDTree==tree)
  north<-m$S.OC;south<-m$J.OC
  val1<-.widths(trw,north,south,height)
  east<-m$V.OC;west<-m$Z.OC
  val2<-.widths(trw,east,west,height)
  range<-c(min(val1$width),min(val2$width),max(val1$width),max(val2$width))
  rangeX<-c((ceiling(max(range)/10)*-10),(ceiling(max(range)/10)*10))
  print(rangeX)
  
  if(method=="Schweingruber"){
    g1<-.plotGraph(plot,tree,height,trw,Schw,withAlometry,direction="North-South",rangeX)
    g2<-.plotGraph(plot,tree,height,trw,Schw,withAlometry,direction="East-West",rangeX)
    lab<-c("Schweingruber - N-S","Schweingruber W-E")
    #plot_grid(g1, g2, labels=lab, ncol = 2, nrow = 1)
  }else{;}
  if(method=="Braam"){
    g1<-.plotGraph(plot,tree,height,trw,Braa,withAlometry,direction="North-South",rangeX)
    g2<-.plotGraph(plot,tree,height,trw,Braa,withAlometry,direction="East-West",rangeX)
    lab<-c("Braam - N-S","Braam W-E")
  }else{;}
  if(method=="Alestalo"){
    g1<-.plotGraph(plot,tree,height,trw,Ales,withAlometry,direction="North-South",rangeX)
    g2<-.plotGraph(plot,tree,height,trw,Ales,withAlometry,direction="East-West",rangeX)
    lab<-c("Alestalo - N-S","Alestalo W-E")
  }else{;}
  plot_grid(g1, g2, labels=lab, ncol = 2, nrow = 1)
}

######################
### Function drawing graph of cross-section profile
###	
### 
######################

drawCrossSectionProfile<-function(trw,plot=1,tree=1,level=1,show.legend=T){
  w <- .seriesTRW_two(trw,tree,plot,level)
  #print(w)
  #w <- .replace(w)
  widths<-.widthsCalculation(w)
  mx<-max(widths[,2:5])
  lim<-c((ceiling(mx/10)*-10),(ceiling(mx/10)*10))
  #lim<-c((ceiling(mx)*-1),ceiling(mx))
  #print(lim)
  widths$S<-widths$S*-1;widths$E<-widths$E*-1
  elipse<-.elipseFun(widths[1,])
  for(i in 2:length(widths$cambAge)){
    elipse2<-.elipseFun(widths[i,])
    elipse<-rbind(elipse,elipse2)
  }
  graph<-ggplot(elipse,aes(x,y)) + geom_point(aes(colour = factor(calYear)),size=0.01)+geom_path(aes(colour = factor(calYear)),size=0.5)+scale_fill_brewer(palette="Spectral")
  graph<-graph + scale_x_continuous(limits = lim, breaks=seq(lim[1],lim[2],10), labels=abs(seq(lim[1],lim[2],10))) + scale_y_continuous(limits = lim, breaks=seq(lim[1],lim[2],10), labels=abs(seq(lim[1],lim[2],10)))
  graph<-graph + labs(colour="Calendar year", x="East     <-  ->     West", y="South     <-  ->     North")
  if(show.legend==T){graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))}
  else{graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")}
  #graph
  return(graph)
}
######################
### Function drawing graph of bai in diferent levels
###	
### 
######################

drawBai<-function(baiFile,plot=1,tree=1,logscale=F,show.legend=T){
  treeForPlot<-subset(.IDdistinct_simple(baiFile),IDPlot==plot & IDTree==tree)
  baiForPlot<-subset(baiFile,select=treeForPlot$Code)
  baiForPlot[baiForPlot == 0] <- NA
  len<-0
  for(i in 1:length(baiForPlot[1,])){
    if(length(na.omit(baiForPlot[,i]))>len)len<-length(na.omit(baiForPlot[,i]))
    #print(c(len,i))
  }
  start<-length(baiForPlot[,1])-len
  end<-length(baiForPlot[,1])
  baiForPlot<-baiForPlot[start:end,]
  a<-baiForPlot[,1]
  
  for(i in 2:length(baiForPlot)){
    a<-c(a,baiForPlot[,i])
  }
  forPlot<-data.frame(x=rep(as.numeric(rownames(baiForPlot)),length(baiForPlot)),
                      y=a,
                      val=rep(1:length(baiForPlot), each=length(as.numeric(rownames(baiForPlot))), length = length(a)))
  if(logscale==F){graph<-ggplot(data=forPlot,aes(x=x, y=y, group = factor(val), colour = factor(val))) + geom_line(size=1.3)}
  else{graph<-ggplot(data=forPlot,aes(x=x, y=y, group = factor(val), colour = factor(val))) + geom_line(size=1.3)+ coord_trans(y = "log10")}
  graph<-graph + labs(colour="Sampling level", x="Calendar year", y="BAI (sq. cm)")
  if(show.legend==T){
    graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }else{
    graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")
  }
  return(graph)
}

######################
### Function drawing graph of taper for entire tree
###	
### 
######################

drawTaper<-function(taperFile,plot=1,tree=1,variant="Taper"){
  #whatDraw<-subset(taperFile,plot==plot & tree==tree)
  whatDraw<-taperFile[taperFile$plot==plot,]
  whatDraw<-whatDraw[whatDraw$tree==tree,]
  if(variant=="Taper"){
    p<-ggplot(data=whatDraw, aes(x=level.height, y=taper, group=1)) + geom_line(size=1.3) + geom_point(size=5)
    p<-p+ labs(colour="Stem taper", x="Sampling height", y="Taper")
  }
  if(variant=="Angle"){
    p<-ggplot(data=whatDraw, aes(x=level.height, y=taper.angle, group=1)) + geom_line(size=1.3) + geom_point(size=5)
    p<-p+ labs(colour="Stem taper", x="Sampling height", y="Taper angle")
  }
  range<-c(0,(ceiling(max(whatDraw$level.height)/10)*10))
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p<-p + scale_x_continuous(limits = range, breaks=whatDraw$level.height, labels=whatDraw$level.height)
  #print(whatDraw)
  return(p)
}

