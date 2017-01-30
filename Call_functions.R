library(dplR) # We need to use dplR package to load RWL data

trw <- read.rwl("f:/Praded_vyska/Tucson/trw_zPASTu.rwl", format="auto") # Test data - GAUK: Apical growth in Praded
meta <- readXL("F:/Praded_vyska/teren.xlsx", rownames=FALSE, header=TRUE, na="", sheet="odbery", stringsAsFactors=TRUE)
po <- readXL("F:/Praded_vyska/teren.xlsx", rownames=FALSE, header=TRUE, na="NA", sheet="p.offset", stringsAsFactors=TRUE)


#Calculations
Excentricity(trw.series=trw)
miss.ring <- EMR(trw, 5, po, method="Both")
apical.series <- apical(trw.series=trw, meta=meta, mr.estimate=miss.ring)



#Graphic functions
drawExcentricityGraph #Doplnit!!
drawCrossSectionProfile(trw=trw,tree=1,plot=1,level=6)
drawAlometryGraph(trw,vyska,plot=1,tree=1,dir="E-W")

