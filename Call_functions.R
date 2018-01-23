library("dplR")
library("ggplot2")
library("cowplot")
library("Rcmdr")
library("xlsx")

#trw <- read.rwl("Last.rwl", format="auto") # Test data - GAUK: Apical growth in Praded
trw <- read.rwl("#Serie_CrossD.rwl", format="auto") # Test data - GAUK: Apical growth in Praded
#meta <- readXL("Last.xlsx", rownames=FALSE, header=TRUE, na="", sheet="odbery",stringsAsFactors=TRUE)
meta <- readXL("Last.xlsx", rownames=FALSE, header=TRUE, na="",sheet="odbery")
po <- readXL("Last.xlsx", rownames=FALSE, header=TRUE, na="", sheet="p.offset",stringsAsFactors=FALSE)

#vysky<-c(97,145,237,293,340,380,531,550,600) #pouze pomocné výšky

#Calculations
exc <- Excentricity(trw.series=trw) # OK
miss.ring <- EMR(trw, 5, po, method="Both") # OK
replacedFile<-replaceMissingRings(trw.series=trw, meta=meta, mr.estimate=miss.ring,no.rings=5,no.series.per.height=4) #OK
apical.series <- apical(trw.series=trw, meta=meta, mr.estimate=miss.ring) # OK
rF<-replacedFile[,1:56]

#taperFile<-taperCalcul(replacedFile,meta) ###!!!### CHYBA - Error in res$taper[i] <- oElipse/treeHeight : replacement has length zero
taperFile<-taperCalcul(rF,meta) # OK 
#baiComplete<-BAIcalculation(replacedFile) # OK, ale výpoèetnì nároèné
baiComplete<-BAIcalculation(rF) #OK

#Graphic functions
drawExcentricityGraph(plot=8,tree=2,heights=meta, trw=rF, exc=excentricity, withAlometry=T, method="Schweingruber") ###OK

drawCrossSectionProfile(replacedFile,plot=8,tree=2,level=2,show.legend=F) #OK

drawBai(baiComplete,plot=8,tree=2,T,T) #OK			

drawTaper(taperFile,plot=8,tree=2,variant="Angle") #OK


