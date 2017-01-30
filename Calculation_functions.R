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

#doplnit o funkci pøipravující data (kdyby mìl uživatel pouze 2 na sebe kolmé vývrty)
BAIcalculation<-function(trw,plot=1,tree=1,level=1){
  w <- .seriesTRW_one(trw,plot=1,tree=1,level=1)
  w <- replace(w)
  widths<-.widthsCalculation(w)
  baiCalculation<-data.frame(date=widths$date,cambAge=widths$cambAge,bai=rep(0,length(widths$date)))
  area<-((pi*abs(widths$W)*abs(widths$N))/4)+((pi*abs(widths$N)*abs(widths$E))/4)+((pi*abs(widths$E)*abs(widths$S))/4)+((pi*abs(widths$S)*abs(widths$W))/4)
  baiCalculation$bai<-area[1]
  for(i in 2:length(baiCalculation$date)){
    j<-i-1
    baiCalculation$bai[i]<-area[i]-area[j]
  }
  
  return(baiCalculation)
}

###########################################
### Function for calculation taper
### Calculates taper for each level
### Used arguments of this function are stored TRW (data.frame name), plot and tree ID (number) and direction of tree core, default is North-South ("N-S"), second Argument is East-West ("E-W")
###########################################

#DODAT funkci na extrakci metadat!!!!!!!!!!!
taperCalculation<-function(trw,heights,plot=1,tree=1,direction="N-S"){
  if(direction=="N-S"){
    serieA<-.extractTRW(trw,plot,tree,"North")
    serieB<-.extractTRW(trw,plot,tree,"South")
  }else{
    serieA<-.extractTRW(trw,plot,tree,"East")
    serieB<-.extractTRW(trw,plot,tree,"West")
  }
  height<-rep(heights, times = length(serie)+1)
  lengths<-.seriesLength(serieA)
  widths<-.calcSumTRW(serieA,lengths)+.calcSumTRW(serieB,lengths)
  curve<-rep(0:(length(serie)), rep(length(serie)+1,length(serie)+1))
  level<-rep(1:(length(serie)+1),length(serie)+1)
  d<-data.frame(curve=curve,level=level,heights=height,widths=widths,taper=rep(0,length(level)))
  taper<-data.frame(first=rep(0,(max(d$level)-1)))
  for(i in 0:(max(d$curve)-1)){
    a<-rep(0,max(d$level)-1)
    b<-subset(d,curve==i)
    for(j in 1:(max(d$level)-1)){
      if(b$widths[j]>0){
        a[j]<-b$widths[j]/(b$heights[max(b$level)]-b$heights[1])
      }else{
        a[j]<-0
      }
    }
    if(i==0){taper <- a}else{taper <- cbind(taper,a)}
  }
  t1<-rep(1,max(d$curve)); t2<-rep("-",(max(d$level-1))); t3<-c(2:max(d$level)); t4 <- paste(t1, t2, t3, sep="")
  rownames(taper)<-t4
  colnames(taper)<-c(1:max(d$curve))
  return(taper)
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
	if (nrow(p.off.2)==1) {value <- round(p.off.2$P.OFFSET/mean, digits=0);
				tab[j,1] <- value;
				tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,]);
				j <- j+1
				} # The subtraction of pith offset an mean tree-ring width is calculated only if there is respective data in pith offset metadata table. Result is appended to output table
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
	if (nrow(p.off.2)==1) {value <- round(pi*(p.off.2$P.OFFSET^2)/mean, digits=0); # The area delimited by the pith offset is divided by mean BAI of last nyrs tree-rings
				tab[j,1] <- value;
				tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,]);
				j <- j+1
				} # The subtraction of pith offset an mean tree-ring width is calculated only if there is respective data in pith offset metadata table. Result is appended to output table
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
	if (nrow(p.off.2)==1) {value <- round(((pi*(p.off.2$P.OFFSET^2)/mean.bai)+(p.off.2$P.OFFSET/mean.trw))/2, digits=0); # BAI and TRW approaches are calculated simultaneously and everaged
				tab[j,1] <- value;
				tab[j,2] <- as.character(data.frame(colnames(trw.series))[i,]);
				j <- j+1
				} # The subtraction of pith offset an mean tree-ring width is calculated only if there is respective data in pith offset metadata table. Result is appended to output table
	 }
}

return(tab)

}

############################################
### Apical growth chronologies
############################################


apical <- function (trw.series, meta, mr.estimate) {
j <- 1
n.ring <- data.frame(ObservedRing=NA, MissingRing=NA, TotalRing=NA, Series=NA)
n.ring.elev <- data.frame(TotalRing=NA, IDPlot=NA, IDTree=NA, IDLevel=NA, Pith=NA)
ID <- .IDdistinct(trw.series)

for (i in (1:ncol(trw.series))){
	df.trw <- data.frame(na.omit(trw.series[,i])) # Getting number of non-NA rings from trw series
	miss.ring <- subset(mr.estimate, subset=Series==as.character(data.frame(colnames(trw))[i,])) # Number of missing rings is taken from dataframe previously created using EMR

	# Appending results to output dataframe
	n.ring[j,1] <- nrow(df.trw) # number of measured rings
	if (nrow(miss.ring)==1) {n.ring[j,2] <- miss.ring[1,1]} # number of missing rings
	n.ring[j,3] <- n.ring[j,2] + n.ring[j,1] # total number of rings (measured+missing)
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
	if (nrow(pith)>0) {ring.estimate <- sum(pith[,"TotalRing"])/nrow(pith);
				n.ring.elev[k,"Pith"] <- "Yes"} 
	if (nrow(pith)==0) {ring.estimate <- sum(subs[,"TotalRing"])/nrow(subs);
				n.ring.elev[k,"Pith"] <- "No"} 

	n.ring.elev[k,"TotalRing"] <- round(ring.estimate, 0) 
	n.ring.elev[k,c("IDPlot", "IDTree", "IDLevel")] <- ((ID$IDAspect)[k,c("IDPlot", "IDTree", "IDLevel")])
	

}

lst <- list(N.ring_Core=n.ring, N.ring_Level=n.ring.elev)
return(lst)
}

