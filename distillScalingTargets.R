

#2k targets from Neukom
cpsEnsemble <- load("/Users/npm4/Dropbox/ChristophTemp12k/DataRaphi/CPS_ensemble.RData")

# Neukom ------------------------------------------------------------------


lats <- seq(-87.5,87.5,by=5)
latbins <- seq(-90,90,by = 30)

for(lb in 1:(length(latbins)-1)){
  fi <- which(lats > latbins[lb] & lats <= latbins[lb+1])
  zonMeans <- apply(ensx[,fi, ,],c(2,3,4),mean)
  lWeights <- cos(lats[fi] * pi/180)
  lWeights <- lWeights/sum(lWeights)
  zmw <- zonMeans * lWeights
  ensMeans <- apply(zmw,c(2,3),sum)
  write.csv(x = cbind(1:2000,t(ensMeans)),file = paste0("~/Dropbox/Temperature12k/",paste0(latbins[lb] ,"to", latbins[lb+1],"ensembleMeansCPS.csv")), row.names = FALSE,col.names = FALSE)
}




#required packages---------------
require(RNetCDF)
require(fields)
library(rgdal)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(parallel)


#reference period
#ref.start<-1911
ref.start <- 1850
#ref.end<-1995
ref.end <- 1900

#start year of the instrumental grid
instr.start<-1850


# HadCrut4 ----------------------------------------------------------------


###
instrfile<-"~/Dropbox/ChristophTemp12k/DataRaphi/input_data/HadCRUT4.3_GraphEM_SP80_18502014_Apr-Mar_corr.nc"
nc<-open.nc(instrfile)
lats<-var.get.nc(nc,1)
lons<-var.get.nc(nc,0)
instr<-var.get.nc(nc,3)
close.nc(nc)
lons.2d<-as.vector(lons)
lats.2d<-as.vector(lats)
lons.grid<-lons
lats.grid<-lats
lons<-unique(lons.2d)
lats<-unique(lats.2d)

nlon<-length(lons)
nlat<-length(lats)
ncells<-length(lons.2d)

#anomalies
instr<-aperm(apply(instr,c(1,2),function(x) x-mean(x[(ref.start-instr.start+1):(ref.end-instr.start+1)])),c(2,3,1))



latbins <- seq(-90,90,by = 30)

for(lb in 1:(length(latbins)-1)){
  fi <- which(lats > latbins[lb] & lats <= latbins[lb+1])
  zonMeans <- apply(instr[,fi ,],c(2,3),mean)
  lWeights <- cos(lats[fi] * pi/180)
  lWeights <- lWeights/sum(lWeights)
  zmw <- zonMeans * lWeights
  ensMeans <- apply(zmw,2,sum)
  out <- cbind(seq(1850,2014),(ensMeans))
  out <- out[-nrow(out),]
  write.csv(x = out,file = paste0("~/Dropbox/Temperature12k/",paste0(latbins[lb] ,"to", latbins[lb+1],"ensembleMeansHadCrut4.csv")), row.names = FALSE,col.names = FALSE)
}




# ERA20C ------------------------------------------------------------------
#required packages---------------
require(RNetCDF)
require(fields)
library(rgdal)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(parallel)

#lat bins
#latbins <- seq(-90,90,by = 30)
#latbins <- c(-90,-30,0,30,90)
#latbins <- c(-90,0,90)
latbins <- c(-90,90)



#reference period
#ref.start<-1911
ref.start <- 1901
#ref.end<-1995
ref.end <- 2000

#start year of the instrumental grid
instr.start<-1900



instrfile<-"~/Dropbox/Temperature12k/t2m_annual_1900_to_2010_era20c.nc"
nc<-open.nc(instrfile)
lats<-var.get.nc(nc,1)
lons<-var.get.nc(nc,2)
years <- var.get.nc(nc,0)
instr<-var.get.nc(nc,3)
close.nc(nc)
lons.2d<-as.vector(lons)
lats.2d<-as.vector(lats)
lons.grid<-lons
lats.grid<-lats
lons<-unique(lons.2d)
lats<-unique(lats.2d)

nlon<-length(lons)
nlat<-length(lats)
ncells<-length(lons.2d)

#anomalies
instr<-aperm(apply(instr,c(1,2),function(x) x-mean(x[(ref.start-instr.start+1):(ref.end-instr.start+1)])),c(2,3,1))





for(lb in 1:(length(latbins)-1)){
  fi <- which(lats > latbins[lb] & lats <= latbins[lb+1])
  zonMeans <- apply(instr[,fi ,],c(2,3),mean)
  lWeights <- cos(lats[fi] * pi/180)
  lWeights <- lWeights/sum(lWeights)
  zmw <- zonMeans * lWeights
  ensMeans <- apply(zmw,2,sum)
  out <- cbind(years,(ensMeans))
  write.csv(x = out,file = paste0("~/Dropbox/Temperature12k/",paste0(latbins[lb] ,"to", latbins[lb+1],"ensembleMeansERA20C.csv")), row.names = FALSE,col.names = FALSE)
}


