###########################################
### Function which separates codes of PlotID, TreeID, Elevation level of corring and aspect (direction of coring) from series name
### Important for future linking of series with metadata + calculation of excentricity indexes for different elevation levels
### The format of series code should be IDPlots_IDTrees_IDElevationAspect
###########################################

.IDdistinct <- function(trw.series, complete=FALSE)
{
  tab <- data.frame(colnames(trw.series), NA, NA, NA, NA) # Creates table with original series names as rows
  colnames(tab) <- c("Original.Code", "IDPlot", "IDTree", "IDLevel", "Aspect") # Appends 4 new columns (IDPlot, ...) to the table
  series.names <- (data.frame(colnames(trw.series)))
  
  for (i in 1:ncol(trw.series)) # Looping one-by-one through the list of series
  { 
    an.series <- as.character(series.names[i,])
    Aspect <- (substr(an.series, nchar(an.series),nchar(an.series))) # Exctracts the last character of original series name (i.e., aspect)...
    tab[i,5] <- Aspect # ... and apends it to the output table.
    split <- strsplit(substr(an.series,1,nchar(an.series)-1), "_") # Splits the rest of the string using _ separator and ...
    unlist <- t(data.frame(as.numeric(unlist(split)))) # ... creates a data frame from it.
    tab[i,2] <- unlist[1,1] # Different parts of data frame are appendend to output table.
    tab[i,3] <- unlist[1,2]
    tab[i,4] <- unlist[1,3]
    rm(an.series, Aspect, split, unlist) # Memory clearing.
  }
  
  S <- subset(tab, subset=Aspect=="S"); colnames(S) <- c("S.OC", "IDPlot", "IDTree", "IDLevel", "S") # Subsetting of table to 4 tables for different aspects
  N <- subset(tab, subset=Aspect=="J"); colnames(N) <- c("J.OC", "J.IDPlot", "J.IDTree", "J.IDLevel", "J")
  V <- subset(tab, subset=Aspect=="V"); colnames(V) <- c("V.OC", "V.IDPlot", "V.IDTree", "V.IDLevel", "V")
  Z <- subset(tab, subset=Aspect=="Z"); colnames(Z) <- c("Z.OC", "Z.IDPlot", "Z.IDTree", "Z.IDLevel", "Z")
  
  if (complete==TRUE)  {
    merge <- merge(S[,1:4],N[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("J.IDPlot", "J.IDTree", "J.IDLevel"), all=T) # Merging different aspects into one file
    merge <- merge(merge,Z[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("Z.IDPlot", "Z.IDTree", "Z.IDLevel"), all=T)
    merge <- merge(merge,V[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("V.IDPlot", "V.IDTree", "V.IDLevel"), all=T) }
  else {
    merge <- merge(S[,1:4],N[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("J.IDPlot", "J.IDTree", "J.IDLevel"), all=F) # Merging different aspects into one file
    merge <- merge(merge,Z[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("Z.IDPlot", "Z.IDTree", "Z.IDLevel"), all=F)
    merge <- merge(merge,V[,1:4], by.x=c("IDPlot", "IDTree", "IDLevel"), by.y=c("V.IDPlot", "V.IDTree", "V.IDLevel"), all=F) }
  
  lst <- list(IDCore=tab, IDAspect=merge)
  return(lst)
}
.IDdistinct_simple<-function(file){
  tab <- data.frame(colnames(file), NA, NA) # Creates table with original series names as rows
  tab<-tab[-c(1),]
  colnames(tab) <- c("Code", "IDPlot", "IDTree") # Appends 4 new columns (IDPlot, ...) to the table
  series.names <- tab$Code
  for (i in 1:length(series.names)) # Looping one-by-one through the list of series
  { 
    an.series <- as.character(series.names[i])
    split <- strsplit(substr(an.series,1,nchar(an.series)-1), "_") # Splits the rest of the string using _ separator and ...
    unlist <- t(data.frame(as.numeric(unlist(split)))) # ... creates a data frame from it.
    tab[i,2] <- unlist[1,1] # Different parts of data frame are appendend to output table.
    tab[i,3] <- unlist[1,2]
    #tab[i,4] <- unlist[1,3]
    rm(an.series, split, unlist) # Memory clearing.
  }
  return(tab)
}#ID distinct for plot and tree only (no aspects etc.)

.extractTreeHeights<-function(meta,plot=1,tree=1){
  tree<-subset(subset(meta,IDPlot==1 & IDTree==1))
  return(tree$Height)
} #extract tree heights from meta data

.extractTRW<-function(trw,plot,tree,direction=""){
  tree<-subset(.IDdistinct(trw)$IDAspect,IDPlot==plot & IDTree==tree)
  ser<-.seriesTRW_one(trw)
  for(i in 2:length(tree$IDLevel)){
    a<-.seriesTRW_one(trw,level=i)
    ser<-cbind(ser,a)
  }
  if(direction==""){
    serie<- ser
  }
  if(direction=="North"){
    serie<- subset(ser, select=tree$S.OC)
  }
  if(direction=="South"){
    serie<- subset(ser, select=tree$J.OC)
  }
  if(direction=="East"){
    serie<- subset(ser, select=tree$V.OC)
  }
  if(direction=="West"){
    serie<- subset(ser, select=tree$Z.OC)
  }
  return(serie)
}# extracts TRW for one chosen tree

######################
### Group of not user-available functions used in excentricity calculation
###	CreateTable to save (i) excentricity indexes and (ii) storing names of input series
###	Caluclation of excentricity indexes according to Schweingruber (1996), Braam et al. (1987) and Alestalo et al. (1971)
######################

.CreateTableTRW<-function(ser,ID){
  table<- data.frame(matrix(nrow=nrow(ser), ncol=(4*nrow(ID)))) # Dataframe to store results ...
  rownames(table) <- rownames(ser) # ... number of rows equal to maximum number of tree-rings in input series
  return(table)
}
.CreateTableMeta<-function(ID){
  table <- data.frame(matrix(nrow=1, ncol=(4*nrow(ID)))) # ... number of columns is 4*number of elevation levels
  return(table)
}

#Functions for calculation of excentricity idexes
.Schwein<-function(ser.s,ser.j,ser.z,ser.v, sub){
  Schwein_ind<- c(sub[,1]/sub[,2], sub[,2]/sub[,1], sub[,3]/sub[,4], sub[,4]/sub[,3]) # Calculates excentrecity index and assignes it to the table
  return(Schwein_ind)
}
.Braam<-function(ser.s,ser.j,ser.z,ser.v, sub){
  # Because we have got 2 perpendicular measurements, I use average of them ! - WE NEED TO FIX THIS FOR SITUATION, WHEN ONLY 3 OR 2 CORES ARE AVAILABLE.# Because we have got 2 perpendicular measurements, I use average of them ! - WE NEED TO FIX THIS FOR SITUATION, WHEN ONLY 3 OR 2 CORES ARE AVAILABLE.
  Braam_ind <- c((sub[,1]-(sub[,3]+sub[,4])/2)/(sub[,1]+(sub[,3]+sub[,4])/2),
                 (sub[,2]-(sub[,3]+sub[,4])/2)/(sub[,2]+(sub[,3]+sub[,4])/2),
                 (sub[,3]-(sub[,1]+sub[,2])/2)/(sub[,3]+(sub[,1]+sub[,2])/2),
                 (sub[,4]-(sub[,1]+sub[,2])/2)/(sub[,4]+(sub[,1]+sub[,2])/2)) 
  return(Braam_ind)
}
.Alestalo<-function(ser.s,ser.j,ser.z,ser.v, sub){
  # Because we have got 2 perpendicular measurements, I use average of them ! - WE NEED TO FIX THIS FOR SITUATION, WHEN ONLY 3 OR 2 CORES ARE AVAILABLE.# Because we have got 2 perpendicular measurements, I use average of them ! - WE NEED TO FIX THIS FOR SITUATION, WHEN ONLY 3 OR 2 CORES ARE AVAILABLE.
  Alest_ind<- c(sub[,1]/(sub[,1]+sub[,2]),
                sub[,2]/(sub[,1]+sub[,2]),
                sub[,3]/(sub[,3]+sub[,4]),
                sub[,4]/(sub[,3]+sub[,4])) 
  return(Alest_ind)
}

######################
### Functions which select names of series or data (excentricity or TRW) for one particular tree
###	Tree selection is made in arguments
### Other functions in this group extracts data from individual data frames or lists
######################

.treeSelect<-function(plot,tree,sourc){
  subset(.IDdistinct(method)$IDAspect,IDPlot==plot & IDTree==tree)
} # Select names of series for particullar tree
.seriesTRW_one<-function(t,plot=1,tree=1,level=1){
  ID<-.IDdistinct(t)
  ID<-subset(ID$IDAspect,IDPlot==plot & IDTree==tree & IDLevel==level)
  #ser.s <- as.character(ID$S.OC) 
  #ser.j <- as.character(ID$J.OC)
  #ser.v <- as.character(ID$V.OC)
  #ser.z <- as.character(ID$Z.OC)
  ser.s <- toString(ID$S.OC[1]) 
  ser.j <- toString(ID$J.OC[1])
  ser.v <- toString(ID$V.OC[1])
  ser.z <- toString(ID$Z.OC[1])
  #print(ser.s)
  w <- t[,c(ser.s, ser.j, ser.v, ser.z)]
  return(w)
}# returns TRW for selected tree including NA values
.seriesTRW_two<-function(t,tree=1,plot=1,level=1){
  #IDa<-.IDdistinct(t)
  #ID<-subset(IDa$IDAspect,IDPlot==plot & IDTree==tree & IDLevel==level)
  #ID<-subset(IDa$IDAspect, IDa$IDAspect$IDPlot==plot & IDa$IDAspect$IDTree==tree & IDa$IDAspect$IDLevel==level)
  #ser.s <- as.numeric(ID["S.OC"])
  #ser.j <- as.numeric(ID["J.OC"])
  #ser.v <- as.numeric(ID["V.OC"])
  #ser.z <- as.numeric(ID["Z.OC"])
  #ser.s <- as.character(ID$S.OC) 
  #ser.j <- as.character(ID$J.OC)
  #ser.v <- as.character(ID$V.OC)
  #ser.z <- as.character(ID$Z.OC)
  #w <- t[,c(ser.s, ser.j, ser.v, ser.z)]
  w <- .seriesTRW_one(t,plot,tree,level)
  l1<-na.omit(w[,1]);l2<-na.omit(w[,2]);l3<-na.omit(w[,3]);l4<-na.omit(w[,4])
  cut<-c((length(w[,1])-length(l1)),(length(w[,1])-length(l2)),(length(w[,1])-length(l3)),(length(w[,1])-length(l4)))
  cut<-min(cut)
  if(cut>0){cut<-cut+1}else{cut<-0}
  w<-w[cut:length(w[,1]),]
  return(w)
}# returns TRW for selected tree without NA values

#Zkontrolovat, jak moc jsou tyto funkce univerzální!!!!
.seriesLength<-function(serie){
  lengths<-c()
  for(i in 1:(length(serie))){
    lengths[i]<-length(na.omit(serie[,i]))
  }
  i<-length(lengths)+1
  lengths[i]<-0
  return(lengths)
}# calculates and return series length
.calcSumTRW<-function(serie,lengths){
  trws<-c()
  apex<-data.frame(apx=rep("NA",length(serie[,1])))
  serie2<-cbind(serie,apex)
  serie3<-cbind(serie,apex)
  for(i in 1:(length(serie2)-1)){
    serie3<-cbind(serie3,serie2)
  }
  
  for(i in 1:length(lengths)){
    x<-length(serie3[,i])
    y<-x-lengths[i]
    if(x==y){
      e<-0
    }else{
      e<-sum(serie3[y:x,i], na.rm = T)
    }
    trws<-c(trws,e)
  }
  return(trws)
}# calculates serie sum for selected series and length 

.widthCalculation<-function(s){
  s[is.na(s)] <- 0
  vypocet<-rep(0,length(s))
  vypocet[1]<-s[1]
  for(i in 2:length(s)){
    vypocet[i]<-vypocet[i-1]+s[i]
  }
  return(vypocet)
} #calculates tree-ring distance from the pith for one serie
.widthsCalculation<-function(width){
  widthsN<-.widthCalculation(width[,1])
  widthsS<-.widthCalculation(width[,2])
  widthsE<-.widthCalculation(width[,3])
  widthsW<-.widthCalculation(width[,4])
  w<-data.frame(date=rownames(width), cambAge=c(1:length(widthsN)),N=widthsN,S=widthsS,E=widthsE,W=widthsW)
  return(w)
}# returns data.frame of distance of each tree ring from the pith for vizualization of - serves mainly for BAI calculation and cross section vizualization
.replace<-function(repl){
  NAN<-0;NAS<-0;NAE<-0;NAW<-0
  ISNAN<-F;ISNAS<-F;ISNAE<-F;ISNAW<-F
  lntghs<-c(na.omit(length(repl[,1])),na.omit(length(repl[,2])),na.omit(length(repl[,3])),na.omit(length(repl[,4])))
  for(i in 1:max(lntghs)){
    if(ISNAN==F){if(is.na(repl[i,1])){NAN<-NAN+1}else{ISNAN<-T}}
    if(ISNAS==F){if(is.na(repl[i,2])){NAS<-NAS+1}else{ISNAS<-T}}
    if(ISNAE==F){if(is.na(repl[i,3])){NAE<-NAE+1}else{ISNAE<-T}}
    if(ISNAW==F){if(is.na(repl[i,4])){NAW<-NAW+1}else{ISNAW<-T}}
  }
  r<-repl
  if(NAN>0){a<-NAN+1;b<-NAN+6;r[1:NAN,1]<-mean(repl[a:b,1])}
  if(NAS>0){a<-NAS+1;b<-NAS+6;r[1:NAS,2]<-mean(repl[a:b,2])}
  if(NAE>0){a<-NAE+1;b<-NAE+6;r[1:NAE,3]<-mean(repl[a:b,3])}
  if(NAW>0){a<-NAW+1;b<-NAW+6;r[1:NAW,4]<-mean(repl[a:b,4])}
  return(r)
}# function for filling missing data - complete tree ring series (with an average of five last known tree rings) - serves mainly for BAI calculation and cross section vizualization

######################
### Set of private functions for creating graph of cross-section profile
###	
### 
######################
.elipseFun <- function(sirky=widths[1,], npoints=100){
  center<-c(0,0)
  t <- seq(0*pi, 2*pi, length.out=npoints)
  q1 <- data.frame(x=rep(0,npoints/4), y=rep(0,npoints/4), calYear=rep(sirky$date,npoints/4), cambAge=rep(sirky$cambAge,npoints/4))
  q2 <- data.frame(x=rep(0,npoints/4), y=rep(0,npoints/4), calYear=rep(sirky$date,npoints/4), cambAge=rep(sirky$cambAge,npoints/4))
  q3 <- data.frame(x=rep(0,npoints/4), y=rep(0,npoints/4), calYear=rep(sirky$date,npoints/4), cambAge=rep(sirky$cambAge,npoints/4))
  q4 <- data.frame(x=rep(0,npoints/4), y=rep(0,npoints/4), calYear=rep(sirky$date,npoints/4), cambAge=rep(sirky$cambAge,npoints/4))
  
  limA<-1;limB<-npoints*0.25
  q1$x <- center[1] + abs(sirky$W)*cos(t[limA:limB])
  q1$y <- center[1] + abs(sirky$N)*sin(t[limA:limB])
  limA<-(npoints*0.25)+1;limB<-npoints*0.5
  q2$x <- center[1] + abs(sirky$E)*cos(t[limA:limB])
  q2$y <- center[1] + abs(sirky$N)*sin(t[limA:limB])
  limA<-(npoints*0.75)+1;limB<-npoints
  q3$x <- center[1] + abs(sirky$W)*cos(t[limA:limB])
  q3$y <- center[1] + abs(sirky$S)*sin(t[limA:limB])
  limA<-(npoints*0.5)+1;limB<-npoints*0.75
  q4$x <- center[1] + abs(sirky$E)*cos(t[limA:limB])
  q4$y <- center[1] + abs(sirky$S)*sin(t[limA:limB])
  
  vystup<-rbind(q1,q2,q4,q3)
  return(vystup)
}# výpočet parametrů elipsy - jednotlivé body

######################
### Set of private functions for creating graph of stem allometry
###	
### 
######################
.seriesCalc<-function(serie,heights){
  lengths<-.seriesLength(serie)
  length2<-c()
  for(i in 1:(length(serie)+1)){
    l<-c(lengths[i:length(lengths)],rep(0,times=i-1))
    length2<-c(length2,l)
  }
  year<-c()
  lab<-as.numeric(rownames(serie))
  yr<-length2[0:(length(serie))+1]
  for(i in length(yr):1){
    a<-length(serie[,1])-yr[i]
    z<-rep(lab[a],times=length(yr))
    year<-c(year,z)
  }
  
  trws<-.calcSumTRW(serie,length2)
  curve<-data.frame(curve=rep(0:(length(serie)), rep(length(serie)+1,length(serie)+1)),
                    year=year,
                    level=rep(1:(length(serie)+1),length(serie)+1),
                    height=rep(heights, times = length(serie)+1),
                    length=length2,
                    trwSum=trws#,
                    #trwAvg=trwa
  )
  return(curve)
}#vrátí data.frame s parametry pro vykreslení alometrického grafu
.dataGraphAlometry<-function(trw,vysky,plot=1,tree=1,dir="N-S"){
  if(dir=="N-S"){
    a<-1
    b<-2
  }else{
    a<-3
    b<-4
  }
  tr<-subset(.IDdistinct(trw)$IDAspect,IDPlot==plot & IDTree==tree)
  ser<-.seriesTRW_one(trw,tree=tree,plot=plot,level=1)
  left<-data.frame(ser[a])
  right<-data.frame(ser[b])
  #print(tree)
  for(i in 2:length(tr$IDLevel)){
    #ser<-.seriesTRW_one(trw,tree=tree,plot=plot,level=i)
    #lvl<-tr$IDLevel[i]
    lvl<-i
    ser<-.seriesTRW_one(trw,tree=tree,plot=plot,level=lvl)
    #print(ser)
    left<-cbind(left,ser[a])
    right<-cbind(right,ser[b])
  }
  l<-.seriesCalc(left,treeHeights);r<-.seriesCalc(right,treeHeights)
  l$trwSum<-l$trwSum*-1;l$length<-l$length*-1
  e<-l$height*-1;e<-c(e,r$heigh)
  LR<-rbind(l,r);LR<-cbind(LR,e)
  return(LR)
}

######################
### Set of private functions for creating graph of excentricity
###	
### 
######################
#set of calculation functions
.widthGraphExcentricity<-function(trw,meta,vysky){
  sirky<-rep(0,length(vysky))
  for(i in 1:length(meta)){
    sirky[i]<-sum(trw[,as.character(meta[i])], na.rm = T)
  }
  sirky[length(vysky)]<-0
  return(data.frame(height=vysky,width=sirky)) 
} #Vypočte vzdálenost od osy y, pro západ a jih (levá cást grafu) se musí vynásobit -1
.dataLine<-function(exct,meta){
  values<-data.frame(exct[,as.character(meta[1])])
  for(i in 2:length(meta)){
    values<-cbind(values,exct[,as.character(meta[i])])
  }
  colnames(values)<-meta
  return(values)
} #Získá data pro excentricitu
.widths<-function(trw,left,right,heights){
  widthsLeft<-.widthGraphExcentricity(trw,left,heights)
  widthsRight<-.widthGraphExcentricity(trw,right,heights)
  widthsLeft$width<-widthsLeft$width*-2
  widthsRight$width<-widthsRight$width*2
  widthsRight<-widthsRight[order(-widthsRight$height),]
  val<-rbind(widthsLeft,widthsRight)
  return(val)
}
.calcSerieRight<-function(serie,vyska){
  serie<-na.omit(serie)
  linie<-data.frame(roky=c(1:length(serie)),hodnoty=serie*10)
  linie$roky<-(linie$roky-1)
  linie$hodnoty<-linie$hodnoty+(vyska-(max(linie$hodnoty)-min(linie$hodnoty))/2)
  return(linie)
} #pro kombinovaný graf vypočte vzdálenosti na ose X a Y pro pravou stranu
.calcSerieLeft<-function(serie,vyska){
  serie<-na.omit(serie)
  linie<-data.frame(roky=c(1:length(serie)),hodnoty=serie*10)
  linie$roky<-(linie$roky-1)*-1
  linie$hodnoty<-linie$hodnoty+(vyska-(max(linie$hodnoty)-min(linie$hodnoty))/2)
  return(linie)
} #pro kombinovaný graf vypočte vzdálenosti na ose X a Y pro levou stranu

#set of drawing functions
.linie<-function(graf, vyska, dataLeft, dataRight){
  line1<-.calcSerieLeft(dataLeft,vyska)
  line2<-.calcSerieRight(dataRight,vyska)
  graf<-graf+geom_line(data=line1, aes_string(line1$roky, line1$hodnoty), size=1, colour="#CC0000")
  graf<-graf+geom_line(data=line2, aes_string(line2$roky, line2$hodnoty), size=1, colour="#CC0000")
  return(graf)
} #vykreslí jednu linii dle parametrů

.plotGraph<-function(plot,tree,heights,trw,exc,withAlometry=T,direction="North-South",rangeX){
  left<-0
  right<-0
  #meta<-.treeSelect(tree,plot,exc)
  met<-subset(.IDdistinct(trw)$IDAspect,IDPlot==plot & IDTree==tree)
  left<-met$S.OC;right<-met$J.OC
  if(direction!="North-South"){left<-met$Z.OC;right<-met$V.OC}
  serieLeft<-.dataLine(exc,left) #extc
  serieRight<-.dataLine(exc,right)
  val<-.widths(trw,left,right,heights)
  #print(val)
  if(withAlometry==T){
    graph<-ggplot(val, aes(width, height)) + geom_point(size = 4, color="#CC0000")+ geom_path(size=1,color="#000099")
    graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    for(i in 1:ncol(serieLeft)){
      graph<-.linie(graph,heights[i],serieLeft[,i],serieRight[,i])
    }
    graph<-graph+scale_x_continuous(name="Stem width (cm)",limits=rangeX, breaks=seq(rangeX[1],rangeX[2],10), labels=(abs(seq(rangeX[1],rangeX[2],10))/10))
    graph<-graph+scale_y_continuous(name="Stem height (cm)",limits=c(0,ceiling(max(heights)/100)*100),breaks=heights,labels=heights)
  }else{
    graph<-ggplot(val, aes(width, height)) + geom_point(size = 4, color="#FFFFFF")+ geom_path(size=1,color="#FFFFFF")
    graph<-graph+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    for(i in 1:ncol(serieLeft)){
      graph<-.linie(graph,heights[i],serieLeft[,i],serieRight[,i])
    }
    graph<-graph+scale_x_continuous(name="Stem width (cm)",limits=rangeX, breaks=seq(rangeX[1],rangeX[2],10), labels=(abs(seq(rangeX[1],rangeX[2],10))/10))
    graph<-graph+scale_y_continuous(name="Stem height (cm)",limits=c(0,ceiling(mx/100)*100),breaks=heights,labels=heights)
  }
  return(graph)
}


.dataLine<-function(exc,met){
  values<-data.frame(exc[as.character(met[1])])
  for(i in 2:length(met)){
    values<-cbind(values,exc[as.character(met[i])])
  }
  colnames(values)<-met
  return(values)
} #Gats data for tree excentricity

