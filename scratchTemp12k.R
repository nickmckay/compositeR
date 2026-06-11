
library(geoChronR) #devtools::install_github("nickmckay/geoChronR")
library(lipdR)
library(purrr)
library(magrittr)
library(ggplot2)
library(compositeR)#devtools::install_github("nickmckay/compositeR")
#load database
D <- readLipd("../lipdFilesWithEnsembles/")

#remove bad values
D$MD97_2151.Assemblage <- NULL
D$NA87_22.Assemblage <- NULL
D$SO42_74KL.Assemblage <- NULL

#extract timeseries
TS <- extractTs(D)


sg <- pullTsVariable(TS,variable = "interpretation1_seasonalityGeneral")
ic <- pullTsVariable(TS,"paleoData_inCompilation")

#filter by compilation and seasonality
tu <- which(tolower(ic) == "temp12kensemble" & (tolower(sg) == "annual" | tolower(sg) == "summeronly" | tolower(sg) == "winteronly"))
fTS <- TS[tu]

ls <- map_dbl(fTS,function(x) sum(!is.na(x$paleoData_values) & !is.na(x$age)))
ls2 <- map_dbl(fTS,function(x) length(x$paleoData_values))

fTS <- fTS[which(ls > 10 & ls2 >10)]


#bin the TS
binvec <-  seq(-50, to = 12050, by = 100)
binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

# #create a matrix of binned data
# binMat <- as.matrix(map_dfc(fTS,simpleBinTs,binvec))
# binMat[is.nan(binMat)] <- NA
#
#
# #create a matrix of binned data pulled from the ensembles
# binMatR <- as.matrix(map_dfc(fTS,sampleEnsembleThenBinTs,binvec))
# binMatR[is.nan(binMatR)] <- NA
#
#
# #let's test it out!
#
# nick <- rerun(50,compositeEnsembles(fTS,binvec,spread = TRUE,duration = 3000, searchRange = c(-50,7000)))
#
# compMatNick <- map_dfc(nick,magrittr::extract)
#
#
#
# cody <- rerun(50,compositeEnsembles(fTS,binvec,stanFun = standardizeOverInterval,interval = c(3000,5000)))
#
# compMatCody <- map_dfc(cody,magrittr::extract)
#
#
# #quick scale to last 2k
# l2k <- which(binAges<2000)
#
# nm2k <- map_dbl(compMatNick[l2k,],mean)
# ns2k <- map_dbl(compMatNick[l2k,],sd)
# cm2k <- map_dbl(compMatCody[l2k,],mean)
# cs2kyyu <- map_dbl(compMatCody[l2k,],sd)
#
# nickScaled <- scale(as.matrix(compMatNick),center = nm2k,scale = ns2k)
# codyScaled <- scale(as.matrix(compMatCody),center = cm2k,scale = cs2k)
#
# #approximate temperature
# nst <- nickScaled/abs(mean(apply(nickScaled,2,diff)[1,])/.5)
# cst <- codyScaled/abs(mean(apply(codyScaled,2,diff)[1,])/.5)
#
# library(ggplot2)
#
# plotTimeseriesEnsRibbons(X = binAges,Y = nst,probs = c(0.25,.5,.75)) %>%
#   plotTimeseriesEnsRibbons(X = binAges,Y = cst,colorHigh = "red",alp = .5,lineColor = "purple",probs = c(0.25,.5,.75))+
#   scale_x_reverse(name = "Age (BP)")+
#   ylab("Approximate temperature scale")
#
#
#
# # Latitudinal gradients ---------------------------------------------------
#
#

#setup ensemble
nens <- 10


#composite lat bins
latbins <- seq(-90,90,by = 30)
lat <- geoChronR::pullTsVariable(fTS,"geo_latitude")

#load in scaling data
targets <- list.files("~/Dropbox/Temperature12k/",pattern = "PAGES2k",full.names = TRUE)
targetsShort <- list.files("~/Dropbox/Temperature12k/",pattern = "PAGES2k",full.names = FALSE)
sw <- 100 #pages 2k scaling window
targ <- purrr::map(targets,read.csv)

scaled <- comps <- counts <- c()
ensOut <- vector(mode = "list",length = nens)
library(foreach)
library(doParallel)
registerDoParallel(4)

allLatBins <- vector(mode = "list",length = 4)
allLatBins[[1]] <- seq(-90,90,by = 30)
allLatBins[[2]] <- c(-90,-30,0,30,90)
allLatBins[[3]] <- c(-90,0,90)
allLatBins[[4]] <- c(-90,90)

#for(alb in 1){
alb <- 1
latbins <- allLatBins[[alb]]
ensOut[[alb]] <- foreach(i = 1:nens) %dopar% {
  scaled <- c()
  for(lb in 1:(length(latbins)-1)){
    fi <- which(lat > latbins[lb] & lat <= latbins[lb+1])
    tc <- compositeEnsembles(fTS[fi],binvec,spread = TRUE,duration = 3000, searchRange = c(0,7000),gaussianizeInput = FALSE)
    #tc <- compositeEnsembles(fTS[fi],binvec,spread = spread,...)
    comps <- cbind(comps,tc$composite)
    counts <- cbind(counts,tc$count)

  #  thisTarget <- which(grepl(targets,pattern = paste0(latbins[lb],"to",latbins[lb+1])))
    thisTarget <- which(stringr::str_starts(string = targetsShort, paste0(latbins[lb] ,"to", latbins[lb+1],"-scaleWindow",sw,"-PAGES2k.csv")))

    if(length(thisTarget) != 1){
      stop("target matching problem")
    }

    thisScaled <- scaleComposite(composite = tc$composite,binvec = binvec,scaleYears = 1950-targ[[thisTarget]][,1],scaleData = targ[[thisTarget]][,-1],scaleWindow = 1950-c(0,2000))


    scaled <- cbind(scaled,thisScaled)
  }

  #weight by areas
  zonalWeights <- sin(latbins[-1]*pi/180)-sin(latbins[-length(latbins)]*pi/180)
  zonalWeights <- zonalWeights/sum(zonalWeights)


  zonalNames <- stringr::str_c(latbins[-1]," to ",latbins[-length(latbins)])
  scaledDf <- as.data.frame(scaled[,-ncol(scaled)])
  names(scaledDf) <- zonalNames
  scaledDf$year <- binAges
  scaledDf$GlobalMean <- rowSums(t(t(scaled)*zonalWeights))
  #scaledDf$counts <- scaled[,ncol(scaled)]



  return(scaledDf)
  # ensOut[[i]] <- scaledDf
  # print(i)
}

save(list = c("ensOut"),file = "~/Dropbox/Temperature12k/12kensOut10.RData")


test <- ensOut[[1]][[1]]
for(i in 2:nens){
  test <- test + ensOut[[1]][[i]]
}
test <- test/nens
test <- test[,c(12,2,4,6,8,10,14)]
names(test) <- c("year",zonalNames)

write_csv(test,"~/Dropbox/Holocene GMST/CPS12k/meanSampleDepth.csv")

#plotting!
library(magrittr)
library(purrr)
allLatMeans <- map_dfc(ensOut[[alb]],extract2,"GlobalMean")
#allCounts <- map_dfc(ensOut[[alb]],extract2,"counts")

#write out data
settings <- paste0(nens,"-",length(latbins)-1,"bands-",sw,"yr2kwindow")
write_csv(path = paste0("~/Dropbox/Holocene GMST/CPS12k/globalMean",settings,".csv"),x = allLatMeans)



for(lb in 1:(length(latbins)-1)){
lbn <- paste0(latbins[lb],"to",latbins[lb+1])
out <- cbind(binAges,as.matrix(map_dfc(ensOut[[alb]],extract2,lb)))
write.csv(file = paste0("~/Dropbox/Holocene GMST/CPS12k/",lbn,settings,".csv"),x = out)

}

library(ggplot2)
globMean <- plotTimeseriesEnsRibbons(X = binAges,Y = as.matrix(allLatMeans),x.bin = seq(-1,12000,by = 10))+
  scale_x_reverse(name = "Year (BP)",breaks = seq(0,12000,2000),oob = scales::squish)+
  scale_y_continuous(name = "Temperature (deg C) (wrt 1000-2000 AD)",limits = c(-10,5),oob = scales::squish)+
  theme_bw()+
  ggtitle("Global Mean Temperature (Composite Plus Scale)")
globMean

ggsave(filename = paste0("~/Dropbox/Temperature12k/GlobalMean12k-2k-",alb,".pdf"),globMean )

#plot bands:


colorsHi <- RColorBrewer::brewer.pal(6,"Dark2")
#colorsHi <- c("blue3", "darkgreen","chocolate4","darkred","darkorchid4","aquamarine4")
#colorsLo <- c("darkslategray1", "darkolivegreen2","chocolate1","coral","darkorchid1","aquamarine")

plot12k <- ggplot()

for(lb in 1:(length(latbins)-1)){

  out <- as.matrix(map_dfc(ensOut[[alb]],extract2,lb))

  plot12k <- plotTimeseriesEnsRibbons(plot12k,X = binAges, Y = out,alp = .5,colorHigh = colorsHi[lb],lineColor = colorsHi[lb],lineWidth = 1,x.bin = seq(-1,12000,by = 10))+
    geom_text(aes(x = 6000), y = (lb * 1.5) - 11 ,label = paste(latbins[lb],"to",latbins[lb+1]),color = colorsHi[lb])

}
plot12k <- plot12k +
  scale_x_reverse(name = "Year (BP)",breaks = seq(0,12000,2000))+
  scale_y_continuous(name = "Temperature (deg C) (wrt 1000-2000 AD)",oob = scales::squish)+
  ggtitle("Zonal Mean Temperature (Composite Plus Scale)")

  theme_bw()

ggsave(filename = paste0("~/Dropbox/Temperature12k/LatBands12k-2k-",alb,".pdf"),plot12k )
plot12k


#
# tidyScale <- tidyr::pivot_longer(scaledDf,cols = -year,names_to = "Latitude Band")
# tidyScale$`Latitude Band` <- factor(tidyScale$`Latitude Band`,levels = rev(zonalNames))
#
# library(ggplot2)
# ggplot()+geom_line(data = tidyScale,aes(x = year, y = value, colour = `Latitude Band`))+
#   geom_line(aes(x = binAges,y = GlobalMean),color = "black",size = 1.5)+
#   scale_x_reverse(name = "Year BP")+
#   scale_y_continuous(name = "Temperature (wrt 0-2000) (deg C)" )+
#   theme_bw()




#plot 2k reconstructions
targets <- list.files("~/Dropbox/Temperature12k/",pattern = "PAGES",full.names = TRUE)
targetsShort <- list.files("~/Dropbox/Temperature12k/",pattern = "PAGES",full.names = FALSE)

targ <- purrr::map(targets,read.csv)

plot2k <- ggplot()


for(lb in 1:(length(latbins)-1)){
plotlb <- ggplot()
  thisTarget <- which(stringr::str_starts(string = targetsShort, paste0(latbins[lb],"to",latbins[lb+1])))

  if(length(thisTarget) != 1){
    stop("target matching problem")
  }

  out <- as.matrix(targ[[thisTarget]])
  out2 <- as.matrix(map_dfc(ensOut[[alb]],extract2,lb))
  ba2 <- 1950-binAges
  out2 <- out2[which(ba2 > 0), ]
  ba2 <- ba2[which(ba2 > 0) ]


  #plot this band
  plotlb <- plotTimeseriesEnsRibbons(plotlb,X = out[,1], Y = scale(out[,-1],scale = FALSE),alp = .8,colorHigh = colorsHi[lb],lineColor = colorsHi[lb],lineWidth = 1,x.bin = seq(0,2000,by=2)) %>%
   #plotTimeseriesEnsRibbons(X = ba2, Y = scale(out2,scale = FALSE),alp = .4,colorHigh = colorsHi[lb],lineColor = colorsHi[lb],lineWidth = 1,x.bin = seq(0,2000,by = 10))+
    plotTimeseriesEnsLines(X = ba2, Y = scale(out2,scale = FALSE))+
    scale_x_continuous(name = "Year (AD)",breaks = seq(0,2000,500))+
    scale_y_continuous(name = "Temperature (deg C) (wrt 1-2000 AD)",oob = scales::squish)+
    ggtitle( paste(latbins[lb],"to",latbins[lb+1]))+
    theme_bw()

  ggsave(filename = paste0("~/Dropbox/Temperature12k/12k2kcompLat_",latbins[lb],"to",latbins[lb+1],"-",alb,".pdf"),plot = plotlb)

  #plot all of them
  plot2k <- plotTimeseriesEnsRibbons(plot2k,X = out[,1], Y = out[,-1],alp = .5,colorHigh = colorsHi[lb],lineColor = colorsHi[lb],lineWidth = 1,x.bin = seq(0,2000,by=2))+
    geom_text(aes(x = 1500), y = (lb * .35), label = paste(latbins[lb],"to",latbins[lb+1]),color = colorsHi[lb])

  out <- as.matrix(map_dfc(ensOut[[alb]],extract2,lb))



}


#
# plot2k+
# scale_x_continuous(name = "Year (AD)",breaks = seq(0,2000,500))+
#   scale_y_continuous(name = "Temperature (deg C) (wrt 1000-2000 AD)",oob = scales::squish)+
#   ggtitle("Zonal Mean Temperature (Composite Plus Scale)")+
# theme_bw()
#
# ggsave(filename = "~/Dropbox/Temperature12k/LatBands2k.pdf",plot2k )

#
# #test scaling
#
# y2k <- targ[[1]]$binAges
# g <- bin(values = as.matrix(targ[[1]][5]),time = as.matrix(y2k),binvec = 1950-seq(-50,950,by=100))
# sd(g$y)
# mean(g$y)
#
# apply(targ[[1]][1901:2000, ],2,sd)
#
#
# apply(test[1:10,],2,sd)

