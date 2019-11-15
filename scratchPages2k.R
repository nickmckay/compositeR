library(geoChronR)
library(lipdR)
library(purrr)
library(magrittr)
#load database
D <- readLipd("~/Dropbox/LiPD/PAGES2k/Temp_v2_1_0/")

#extract timeseries
TS <- extractTs(D)

#filter by compilation
fTS <- filterTs(TS,"paleoData_useInGlobalTemperatureAnalysis ==  TRUE")

ls <- map_dbl(fTS,function(x) sum(!is.na(x$paleoData_values) & !is.na(x$year)))
ls2 <- map_dbl(fTS,function(x) length(x$paleoData_values))

fTS <- fTS[which(ls > 10 & ls2 >10)]


#bin the TS
binvec <-  seq(0, to = 2000, by = 1)
binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))


#composite lat bins
latbins <- seq(-90,90,by = 30)
lat <- geoChronR::pullTsVariable(fTS,"geo_latitude")

#load in scaling data
targets <- list.files("~/Dropbox/Temperature12k/",pattern = "ERA",full.names = TRUE)
targ <- purrr::map(targets,read.csv)


library(foreach)
library(doParallel)
registerDoParallel(6)
#set up ensembles?
nens <- 50

scaled <- comps <- counts <- c()
foreach(lb = 1:(length(latbins)-1)) %dopar% {
  scaleEns <- c()
  fi <- which(lat > latbins[lb] & lat <= latbins[lb+1])
  for(n in 1:nens){
  tc <- compositeEnsembles(fTS[fi],binvec,ageVar = "year",spread = TRUE,duration = 100, searchRange = c(0,2000),binFun = simpleBinTs)
  comps <- cbind(comps,tc$composite)
  counts <- cbind(counts,tc$count)

  thisTarget <- which(grepl(targets,pattern = paste0(latbins[lb],"to",latbins[lb+1])))
  if(length(thisTarget) != 1){
    stop("target matching problem")
  }

  thisScaled <- scaleComposite(composite = tc$composite,binvec = binvec,scaleYears = targ[[thisTarget]][,1],scaleData = targ[[thisTarget]][,-1])


  scaled <- cbind(scaled,thisScaled)
  scaleEns <- cbind(scaleEns,thisScaled)
  }
  out <- cbind(binAges,scaleEns)
  write.csv(x = out,file = paste0("~/Dropbox/Temperature12k/",paste0(latbins[lb] ,"to", latbins[lb+1],"PAGES2k.csv")), row.names = FALSE,col.names = FALSE)
}



#plot 2k reconstructions
targets <- list.files("~/Dropbox/Temperature12k/",pattern = "PAGES",full.names = TRUE)
targ <- purrr::map(targets,read.csv)

colorsHi <- RColorBrewer::brewer.pal(6,"Set3")
#colorsHi <- c("blue3", "darkgreen","chocolate4","darkred","darkorchid4","aquamarine4")
#colorsLo <- c("darkslategray1", "darkolivegreen2","chocolate1","coral","darkorchid1","aquamarine")

plot2k <- ggplot()

for(lb in 1:(length(latbins)-1)){

  thisTarget <- which(grepl(targets,pattern = paste0(latbins[lb],"to",latbins[lb+1])))
  if(length(thisTarget) != 1){
    stop("target matching problem")
  }

  out <- as.matrix(targ[[thisTarget]])

  plot2k <- plotTimeseriesEnsRibbons(plot2k,X = out[,1], Y = out[,-1],alp = .5,colorHigh = colorsHi[lb],lineColor = colorsHi[lb],lineWidth = 1)
}
plot2k





scaleEns <- as.matrix(out[,-1])
plotTimeseriesEnsLines(X = binAges,Y = scaleEns,maxPlotN = 100)





#weight by areas
zonalWeights <- sin(latbins[-1]*pi/180)-sin(latbins[-length(latbins)]*pi/180)
zonalWeights <- zonalWeights/sum(zonalWeights)


zonalNames <- stringr::str_c(latbins[-1]," to ",latbins[-length(latbins)])
scaledDf <- as.data.frame(scaled)
names(scaledDf) <- zonalNames
scaledDf$year <- binAges
GlobalMean <- rowSums(t(t(scaled)*zonalWeights))

tidyScale <- tidyr::pivot_longer(scaledDf,cols = -year,names_to = "Latitude Band")
tidyScale$`Latitude Band` <- factor(tidyScale$`Latitude Band`,levels = rev(zonalNames))

library(ggplot2)
ggplot()+geom_line(data = tidyScale,aes(x = year, y = value, colour = `Latitude Band`))+
  geom_line(aes(x = binAges,y = GlobalMean),color = "black",size = 1.5)+
  scale_x_continuous(name = "Year BP")+
  scale_y_continuous(name = "Temperature (wrt 1850-2000) (deg C)")+
  theme_bw()




