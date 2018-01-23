###########################################
### Function for calculation of eccentricity index according to Schweingruber (1996), Braam et al. (1987) and Alestalo et al. (1971)
### For details about eqations used for calculation of indexes, see notes below or summary in Tumajer, Treml (2013), Geochronometria 40(1): 59-76.
### Procedure: [1] creates new output tabs (private functions) [2] splits table with RWL series into 4 tables according aspect (data taken from .IDformat$Aspect - private function) [3] calculates indexes and stores them in three separate tables -> list
### Stores a results (i.e., eccentricity indexes) in RWL (in R) format - other functions of eg., dplR can be used in further calculations
### Title of the eccentricity series name: IDPlot_IDTree_IDLevelOrientation, where Orientation means Direction from which eccentricity index is calculated
###########################################

Excentricity <- function (trw.series, complete=FALSE)
{
  ID<-.IDdistinct(trw.series, complete)
  
  tab.exc.schw <-.CreateTableTRW(trw.series,ID$IDAspect)
  colnames.store.schw <- .CreateTableMeta(ID$IDAspect)
  
  tab.exc.braam <-.CreateTableTRW(trw.series,ID$IDAspect)
  colnames.store.braam <- .CreateTableMeta(ID$IDAspect)
  
  tab.exc.alestalo <-.CreateTableTRW(trw.series,ID$IDAspect)
  colnames.store.alestalo <- .CreateTableMeta(ID$IDAspect)
  ###################################################################################################
  
  for (i in 1:nrow(ID$IDAspect))
  {
    
    ser.s <- as.character(ID$IDAspect[i,"S.OC"]) # Extracts names of series from the same height level from IDformat$ID Aspect
    ser.j <- as.character(ID$IDAspect[i,"J.OC"])
    ser.z <- as.character(ID$IDAspect[i,"Z.OC"])
    ser.v <- as.character(ID$IDAspect[i,"V.OC"])
    
    sub <- trw.series[,c(ser.s, ser.j, ser.z, ser.v)] # Selects series of those names
    
    tab.exc.schw[,c(4*i-3, 4*i-2, 4*i-1, 4*i)]<-.Schwein(ser.s,ser.j,ser.z,ser.v, sub)
    tab.exc.braam[,c(4*i-3, 4*i-2, 4*i-1, 4*i)]<-.Braam(ser.s,ser.j,ser.z,ser.v, sub)
    tab.exc.alestalo[,c(4*i-3, 4*i-2, 4*i-1, 4*i)]<-.Alestalo(ser.s,ser.j,ser.z,ser.v, sub)
    
    colnames.store.schw[,c(4*i-3, 4*i-2, 4*i-1, 4*i)] <- c(ser.s, ser.j, ser.z, ser.v)
    colnames.store.braam[,c(4*i-3, 4*i-2, 4*i-1, 4*i)] <- c(ser.s, ser.j, ser.z, ser.v)
    colnames.store.alestalo[,c(4*i-3, 4*i-2, 4*i-1, 4*i)] <- c(ser.s, ser.j, ser.z, ser.v)
    
    rm(ser.v, ser.j, ser.z, ser.s) # Clearing memory
  }
  colnames(tab.exc.schw) <- colnames.store.schw[1,]
  colnames(tab.exc.braam) <- colnames.store.braam[1,]
  colnames(tab.exc.alestalo) <- colnames.store.alestalo[1,]
  
  return(list(Schweingruber=tab.exc.schw, Braam=tab.exc.braam, Alestalo=tab.exc.alestalo))
}

###########################################
### Function for calculation of BAI index according to ABC (XYZ),
### A new approach for BAI calculation baseo on approximation of tree ring to elipse
### Used arguments of this function are stored TRW (data.frame name), plot and tree ID (number)
###########################################

#doplnit o funkci pøipravující data (kdyby mìl uživatel pouze 2 na sebe kolmé vývrty??)
BAIcalculation<-function(trw.series){
  IDs<-.IDdistinct(trw.series)
  result<-data.frame(cambAge=c(1:length(rownames(trw.series))))
  rownames(result)<-rownames(trw.series)
  coln<-"cambAge"
  for(i in 1:length(IDs$IDAspect$IDPlot)){
    plot<-IDs$IDAspect$IDPlot[i]
    tree<-IDs$IDAspect$IDTree[i]
    level<-IDs$IDAspect$IDLevel[i]
    #w <- .seriesTRW_two(trw,plot,tree,level)
    w <- .seriesTRW_one(trw.series,plot,tree,level)
    if(length(w[,1])!=length(w[,2]) | length(w[,1])!=length(w[,3]) | length(w[,1])!=length(w[,4])){w <- .replace(w)}
    widths<-.widthsCalculation(w)
    baiCalculation<-data.frame(date=widths$date,cambAge=widths$cambAge,bai=rep(0,length(widths$date)))
    area<-((pi*abs(widths$W)*abs(widths$N))/4)+((pi*abs(widths$N)*abs(widths$E))/4)+((pi*abs(widths$E)*abs(widths$S))/4)+((pi*abs(widths$S)*abs(widths$W))/4)
    baiCalculation$bai<-area[1]
    for(i in 2:length(baiCalculation$date)){
      j<-i-1
      baiCalculation$bai[i]<-area[i]-area[j]
    }
    coln<-c(coln,paste(plot, tree, level, sep="_"))
    b<-c(rep("<NA>",length(result[,1])-length(baiCalculation[,3])),baiCalculation[,3])
    result<-cbind(result, b)
  }
  colnames(result)<-coln
  result$cambAge <- NULL
  return(result)
}

###########################################
### Function for calculation taper
### Calculates taper for each level
### Used arguments of this function are stored TRW (data.frame name), plot and tree ID (number) and direction of tree core, default is North-South ("N-S"), second Argument is East-West ("E-W")
###########################################

taperCalcul<-function(trw.series,meta){
  
  IDa<-.IDdistinct(trw.series)
  ID<-IDa$IDAspect
  res<-data.frame(plot=ID$IDPlot,tree=ID$IDTree,level=ID$IDLevel,level.height=rep(0,length(ID$IDPlot)),taper=rep(0,length(ID$IDPlot)),taper.angle=rep(0,length(ID$IDPlot)))
  
  for(i in 1:length(ID$IDPlot)){
    IDB<-ID[i,]
    j<-IDB$IDLevel+1
    subMeta1<-subset(meta,Plocha_ID==IDB$IDPlot & Strom_ID==IDB$IDTree & Vyska_ID==IDB$IDLevel)
    sb<-length(subset(meta,Plocha_ID==IDB$IDPlot & Strom_ID==IDB$IDTree))+1
    if(j>sb)j<-999
    subMeta2<-subset(meta,Plocha_ID==IDB$IDPlot & Strom_ID==IDB$IDTree & Vyska_ID==j)
    treeHeight1<-subMeta1$Vyska_cm
    treeHeight2<-subMeta2$Vyska_cm
    treeHeight<-as.numeric(treeHeight2)-as.numeric(treeHeight1)
    
    ser.s <- toString(IDB$S.OC)
    ser.j <- toString(IDB$J.OC)
    ser.v <- toString(IDB$V.OC)
    ser.z <- toString(IDB$Z.OC)
    w <- trw.series[,c(ser.s, ser.j, ser.v, ser.z)]
    if(length(w[,1])!=length(w[,2]) | length(w[,1])!=length(w[,3]) | length(w[,1])!=length(w[,4])){w <- .replace(w)}
    widths<-.widthsCalculation(w)
    last<-length(widths$cambAge)
    oNW<-(pi*sqrt(2*(widths$N[last]*widths$N[last]+widths$W[last]*widths$W[last])))/4
    oSW<-(pi*sqrt(2*(widths$S[last]*widths$S[last]+widths$W[last]*widths$W[last])))/4
    oNE<-(pi*sqrt(2*(widths$N[last]*widths$N[last]+widths$E[last]*widths$E[last])))/4
    oSE<-(pi*sqrt(2*(widths$S[last]*widths$S[last]+widths$E[last]*widths$E[last])))/4
    oElipse1<-oNW+oSW+oNE+oSE
    
    if(j!=999){
      IDC<-ID[j,]
      ser.s <- toString(IDC$S.OC)
      ser.j <- toString(IDC$J.OC)
      ser.v <- toString(IDC$V.OC)
      ser.z <- toString(IDC$Z.OC)
      w <- trw.series[,c(ser.s, ser.j, ser.v, ser.z)]
      if(length(w[,1])!=length(w[,2]) | length(w[,1])!=length(w[,3]) | length(w[,1])!=length(w[,4])){w <- .replace(w)}
      widths<-.widthsCalculation(w)
      last<-length(widths$cambAge)
      oNW<-(pi*sqrt(2*(widths$N[last]*widths$N[last]+widths$W[last]*widths$W[last])))/4
      oSW<-(pi*sqrt(2*(widths$S[last]*widths$S[last]+widths$W[last]*widths$W[last])))/4
      oNE<-(pi*sqrt(2*(widths$N[last]*widths$N[last]+widths$E[last]*widths$E[last])))/4
      oSE<-(pi*sqrt(2*(widths$S[last]*widths$S[last]+widths$E[last]*widths$E[last])))/4
      oElipse2<-oNW+oSW+oNE+oSE
    }else{
      oElipse2<-0
    }
    
    
    #print(subMeta)
    #print(c(sb,i,j,treeHeight1,treeHeight2,oElipse1,oElipse2))
    res$level.height[i]<-treeHeight1
    res$taper[i]<-(oElipse1-oElipse2)/treeHeight
    res$taper.angle[i]<-(atan(0.5*res$taper[i])* 180) / (pi)
    #print(res[i,])
  }
  
  return(res)
}

###########################################
### Function for estimating the number of missing rings from pith offset (as metadata tabel)
### There are three options of algorithm for estimating number of missing rings (according to Altman et al. (2016): Forest Ecology and Management, 380:82-89)
### 1] pith offset is divided by the average TRW of last nyrs tree-rings
### 2] pith offset is converted to area (pi*PO^2) and then divided by average BAI of last nyrs tree-rings 
### 3] mean of 1] and 2] (adviced by Altman et al. (2016): Forest Ecology and Management, 380:82-89 as the aproach with the lowest bias)
###########################################

EMR <- function(trw.series, nyrs, p.off, method="TRW"){
  j<-1
  tab <- data.frame(MissingRings=NA, Series=NA)
  
  ### TRW based estimate of number of missing rings
  #########
  if (method=="TRW") {
    for (i in (1:ncol(trw.series))){
      df <- data.frame(na.omit(trw.series[,i])) # trw series with removed NAs
      mean <- mean(df[(1:nyrs),1]) # mean of last nyrs tree-rings 
      
      p.off.2 <- subset(p.off, subset=ID==as.character(data.frame(colnames(trw.series))[i,])) # Subset of pit.offset table to contain only row with tree, which is currently analysed
      
      if (as.numeric(p.off.2["P.OFFSET"])==0 && !is.na(p.off.2["P.OFFSET"])) {tab[j,1] <- 0; tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,])} # If core goes through pith -> 0 missing rings
      else
      {if (!is.na(p.off.2["P.OFFSET"])) {value <- round(as.numeric(p.off.2["P.OFFSET"])/mean, digits=0); # If core misses the pith -> rounded(offset/meanTRW) missing rings
      tab[j,1] <- value;
      }} 
    
    tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,]);
    j <- j+1
    }
  }
  
  ### BAI based estimate of number of missing rings
  #########
  if (method=="BAI") {
    bai <- bai.in(trw.series, subset(p.off, select=c(ID,P.OFFSET))) # Calculation of BAI
    
    # Almost the same script as in case of TRW
    for (i in (1:ncol(bai))){
      df <- data.frame(na.omit(bai[,i])) # bai series with removed NAs
      mean <- mean(df[(1:nyrs),1]) # mean of last nyrs tree-rings 
      
      p.off.2 <- subset(p.off, subset=ID==as.character(data.frame(colnames(bai))[i,])) # Subset of pit.offset table to contain only row with tree, which is currently analysed
      
      if (as.numeric(p.off.2["P.OFFSET"])==0 && !is.na(p.off.2["P.OFFSET"])) {tab[j,1] <- 0; tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,])}
      else
      {if (!is.na(p.off.2["P.OFFSET"])) {value <- round(pi*(as.numeric(p.off.2["P.OFFSET"])^2)/mean, digits=0); # The area delimited by the pith offset is divided by mean BAI of last nyrs tree-rings
      tab[j,1] <- value;
      }} 
      
    tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,]);
    j <- j+1
    }
  }
  
  ### Altman et al. (2016) approach - mean of missing rings estimated based on TRW and BAI
  #########
  if (method=="Both") {
    bai <- bai.in(trw.series, subset(p.off, select=c(ID,P.OFFSET))) # Calculation of BAI
    
    for (i in (1:ncol(bai))){
      df.trw <- data.frame(na.omit(trw.series[,i])) # trw and BAI series with removed NAs
      df.bai <- data.frame(na.omit(bai[,i])) 
      mean.trw <- mean(df.trw[(1:nyrs),1]) # mean of last nyrs tree-rings 
      mean.bai <- mean(df.bai[(1:nyrs),1])
      
      p.off.2 <- subset(p.off, subset=ID==as.character(data.frame(colnames(trw))[i,])) # Subset of pit.offset table to contain only row with tree, which is currently analysed
      
      if (as.numeric(p.off.2["P.OFFSET"])==0 && !is.na(p.off.2["P.OFFSET"])) {tab[j,1] <- 0; tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,])}
      else
      {if (!is.na(p.off.2["P.OFFSET"])) {value <- round(((pi*(as.numeric(p.off.2["P.OFFSET"])^2)/mean.bai)+(as.numeric(p.off.2["P.OFFSET"])/mean.trw))/2, digits=0); # BAI and TRW approaches are calculated simultaneously and everaged
      tab[j,1] <- value;
      }} 
      
     tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,]);
     j <- j+1
    }
  }
  
  return(tab)
  
}

###########################################
### Function for replacing missing tree rings
### Several options is necessearry to fill
### trw.series - file that contains TRW series
### nyrs - number of years fo interpolation
### p.off - pith offset - an estimated distance to the pith
### method - interpolation method
### no.series.per.height - number of cores from the same tree height for interpolation
###########################################
RMR <- function (trw.series, meta, mr.estimate, nyrs=5, nsph=4) {
  j <- 1
  n.ring <- data.frame(ObservedRing=NA, MissingRing=NA, TotalRing=NA, Series=NA, IDCore=NA)
  n.ring.elev <- data.frame(TotalRing=NA, IDCore=NA, Pith=NA)
  ID <- .IDdistinct(trw.series)
  trw.series.descend<-trw.series[ order(-as.numeric(row.names(trw.series))), ]
  
  #One core per sampling height
  if (nsph==1){
    mr.estimate[is.na(mr.estimate)]<-0
    trw.result<-data.frame(matrix(NA, ncol = length(trw.series), nrow = (length(trw.series)+max(mr.estimate$MissingRings))))
    names(trw.result)<-names(trw.series)
    for (i in 1:length(trw.series)){
      ser<-na.omit(trw.series[,i])
      avg<-mean(ser[1:nyrs])
      miss.i<-which(mr.estimate$Series==names(trw.series[i]))
      miss<-mr.estimate$MissingRings[miss.i]           
      ser2<-c(na.omit(trw.series.descend[,i]),rep(avg,miss))
      trw.result[1:length(ser2),i]<-ser2
    }
  }
  
  #Two cores per sampling height
  if(nsph==2){
    mr.estimate[is.na(mr.estimate)]<-0
    n <- substring(names(trw.series),1,(nchar(names(trw.series))-1))
    mr.estimate$IDCore<-n
    n <- unique(n)
    len<-NULL
    for(i in 1:length(trw.series)){
      len<-c(len,length(na.omit(trw.series[,i])))
    }
    mr.estimate$Length<-len
    #sem se musí dodat info o pith presence
    for(i in 1:length(n)){
      ind<-which(n[i]==mr.estimate$IDCore)
      miss<-mr.estimate$MissingRings[ind]
      len<-mr.estimate$Length[ind]
      m<-mean(c((miss[1]+len[1]),(miss[2]+len[2])))
      mr.estimate$LengthTotal[ind]<-m
    }
    trw.result<-data.frame(matrix(NA, ncol = length(trw.series), nrow = max(mr.estimate$LengthTotal)))
    names(trw.result)<-names(trw.series)
    for(i in 1:length(n)){
      ind<-which(n[i]==mr.estimate$IDCore)
      sA<-mr.estimate$Series[ind[1]]
      sB<-mr.estimate$Series[ind[2]]
      serA<-na.omit(trw.series[,sA])
      serB<-na.omit(trw.series[,sB])
      avgA<-mean(serA[1:nyrs])
      avgB<-mean(serB[1:nyrs])
      mA<-mr.estimate$LengthTotal[ind[1]]-mr.estimate$Length[ind[1]]
      mB<-mr.estimate$LengthTotal[ind[2]]-mr.estimate$Length[ind[2]]
      serA2<-c(na.omit(trw.series.descend[,sA]),rep(avgA,mA))
      serB2<-c(na.omit(trw.series.descend[,sB]),rep(avgB,mB))
      trw.result[1:length(serA2),sA]<-serA2
      trw.result[1:length(serB2),sB]<-serB2
    }
  }
  
  #Four cores per sampling height
  if(nsph==4){
    for (i in (1:ncol(trw.series))){
      df.trw <- data.frame(na.omit(trw.series[,i])) # Getting number of non-NA rings from trw series
      miss.ring <- subset(mr.estimate, subset=Series==as.character(data.frame(colnames(trw))[i,])) # Number of missing rings is taken from dataframe previously created using EMR
      
      # Appending results to output dataframe
      n.ring[j,1] <- nrow(df.trw) # number of measured rings
      if (nrow(miss.ring)==1) {n.ring[j,2] <- miss.ring[1,1]} # number of missing rings
      if (!is.na(miss.ring[1,1])) {n.ring[j,3] <- n.ring[j,2] + n.ring[j,1]} else {n.ring[j,3] <- n.ring[j,1]}  # total number of rings (measured+missing)
      n.ring[j,4] <- as.character(data.frame(colnames(trw))[i,]) # name of series
      j <- j+1
    }
    n.ring$IDCore <- substring(n.ring$Series,1,(nchar(n.ring$Series)-1))
    n.ring$MissingRing[is.na(n.ring$MissingRing)]<-0
    for (k in (1:nrow(ID$IDAspect))) {
      ser.s <- as.character((ID$IDAspect)[k,"S.OC"])
      ser.j <- as.character((ID$IDAspect)[k,"J.OC"])
      ser.v <- as.character((ID$IDAspect)[k,"V.OC"])
      ser.z <- as.character((ID$IDAspect)[k,"Z.OC"])
      
      subs <- subset(n.ring, subset=(Series==ser.s|Series==ser.j|Series==ser.v|Series==ser.z)) # subset of estimated number of tree-rings for cores from the same level 
      pith <- subset(subs, subset=MissingRing==0) # subset of number of tree-rings from cores from the same level, which hit the pit
      ### If at least one core goes through the pit, then it is used. Othervise, we use mean of estimated number of rings based on all available cores.
      if (nrow(pith)>0) {ring.estimate <- sum(pith[,"TotalRing"])/nrow(pith);
      n.ring.elev[k,"Pith"] <- "Yes"} 
      if (nrow(pith)==0) {ring.estimate <- sum(subs[,"TotalRing"])/nrow(subs);
      n.ring.elev[k,"Pith"] <- "No"} 
      
      n.ring.elev[k,"TotalRing"] <- round(ring.estimate, 0) 
      cr<-c(paste(((ID$IDAspect)[k,c("IDPlot", "IDTree", "IDLevel")]), collapse="_"))
      #print(cr)
      n.ring.elev[k,c("IDCore")] <- cr
    }
    for (i in 1:length(n.ring.elev$TotalRing)){
      ind<-which(as.character(n.ring$IDCore)==as.character(n.ring.elev$IDCore[i]))
      n.ring$TotalRing[ind]<-n.ring.elev$TotalRing[i]
    }
    n.ring$MissingRing<-n.ring$TotalRing-n.ring$ObservedRing
    n.ring$MissingRing[n.ring$MissingRing<0]<-0
    trw.result<-data.frame(matrix(NA, ncol = length(trw.series), nrow = max(n.ring$TotalRing)))
    names(trw.result)<-names(trw.series)
    for (i in 1:length(trw.series)){
      ser<-na.omit(trw.series[,i])
      avg<-mean(ser[1:nyrs])
      miss.i<-which(n.ring$Series==names(trw.series[i]))
      miss<-n.ring$MissingRing[miss.i]           
      ser2<-c(na.omit(trw.series.descend[,i]),rep(avg,miss))
      trw.result[1:length(ser2),i]<-ser2
    }
  }
  
  start<-as.numeric(rownames(trw.series.descend)[1])
  end<-as.numeric(rownames(trw.series.descend)[1])-length(trw.result[,1])+1
  years<-c(start:end)
  rownames(trw.result)<-years
  trw.result<-trw.result[ order(as.numeric(row.names(trw.result))), ]
  #return(n.ring)
  #return(n.ring.elev)
  return(trw.result)
}

############################################
### Apical growth chronologies
############################################


apical <- function (trw.series, meta, mr.estimate) {
  j <- 1
  n.ring <- data.frame(ObservedRing=NA, MissingRing=NA, TotalRing=NA, Series=NA)
  n.ring.elev <- data.frame(TotalRing=NA, IDPlot=NA, IDTree=NA, IDLevel=NA, Pith=NA, Height.cm=NA)
  ID <- .IDdistinct(trw.series)
  
  for (i in (1:ncol(trw.series))){
    df.trw <- data.frame(na.omit(trw.series[,i])) # Getting number of non-NA rings from trw series
    miss.ring <- subset(mr.estimate, subset=Series==as.character(data.frame(colnames(trw))[i,])) # Number of missing rings is taken from dataframe previously created using EMR
    
    # Appending results to output dataframe
    n.ring[j,1] <- nrow(df.trw) # number of measured rings
    if (nrow(miss.ring)==1) {n.ring[j,2] <- miss.ring[1,1]} # number of missing rings
    if (!is.na(miss.ring[1,1])) {n.ring[j,3] <- n.ring[j,2] + n.ring[j,1]} else {n.ring[j,3] <- n.ring[j,1]}  # total number of rings (measured+missing)
    n.ring[j,4] <- as.character(data.frame(colnames(trw))[i,]) # name of series
    j <- j+1
  }
  
  for (k in (1:nrow(ID$IDAspect))) {
    ser.s <- as.character((ID$IDAspect)[k,"S.OC"])
    ser.j <- as.character((ID$IDAspect)[k,"J.OC"])
    ser.v <- as.character((ID$IDAspect)[k,"V.OC"])
    ser.z <- as.character((ID$IDAspect)[k,"Z.OC"])
    
    subs <- subset(n.ring, subset=(Series==ser.s|Series==ser.j|Series==ser.v|Series==ser.z)) # subset of estimated number of tree-rings for cores from the same level 
    pith <- subset(subs, subset=MissingRing==0) # subset of number of tree-rings from cores from the same level, which hit the pit
    ### If at least one core goes through the pit, then it is used. Othervise, we use mean of estimated number of rings based on all available cores.
    if (nrow(pith)>0) {ring.estimate <- median(pith[,"TotalRing"]);
    n.ring.elev[k,"Pith"] <- "Yes"} 
    if (nrow(pith)==0) {ring.estimate <- median(subs[,"TotalRing"]);
    n.ring.elev[k,"Pith"] <- "No"} 
    
    n.ring.elev[k,"TotalRing"] <- round(ring.estimate, 0) 
    n.ring.elev[k,c("IDPlot", "IDTree", "IDLevel")] <- ((ID$IDAspect)[k,c("IDPlot", "IDTree", "IDLevel")])
    
    meta.subs <- subset(meta, subset=(meta[,1]==n.ring.elev[k,"IDPlots"] & meta[,2]==n.ring.elev[k,"IDTree"] & meta[,3]==n.ring.elev[k,"IDLevel"])
    if (nrow(meta.subs)==1) {n.ring.elev[k,"Height.cm"] <- meta.subs[1,4]}
    
  }
  
  lst <- list(N.ring_Core=n.ring, N.ring_Level=n.ring.elev)
  return(lst)
}


