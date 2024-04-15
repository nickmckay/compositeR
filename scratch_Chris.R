devtools::load_all()
library(lipdR)
library(magrittr)
library(geoChronR)
library(tidyverse)
library(sf)
library(dggridR)
# function to create "square" regions
# defineRegion <- function(latmin,latmax,lonmin,lonmax){
#   nlon = length(seq(lonmin,lonmax,1))
#   nlat = length(seq(latmin,latmax,1))
#   Poly_Coord_df <- data.frame(lon= c(seq(lonmin,lonmax,1),rev(seq(lonmin,lonmax,1))),
#                               lat= c(rep(latmin,nlon),rep(latmax,nlon))
#   ) %>%
#     sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
#     summarise((geometry = sf::st_combine(geometry))) %>%
#     sf::st_cast("POLYGON")
#   return(Poly_Coord_df)
# }
#set arctic circle to plot as panel boundary
arcticCircle <- data.frame(id="A",lon=seq(-180,180), lat= rep(58,length(seq(-180,180))))%>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise((geometry = sf::st_combine(geometry))) %>%
  st_cast("MULTILINESTRING")
#get arctic polygons
arctic <- defineRegion(58,90,-180,180)
# polygonsArctic  <- sf::st_intersection(polygons,arctic)  %>%
#   mutate(lon=st_coordinates(st_centroid(.))[,1]) %>%
#   mutate(lat=st_coordinates(st_centroid(.))[,2]) %>%
#   mutate(area_km2 = units::drop_units(st_area(.)/1000000)) #km^2
#countries to plot
countries <- st_as_sf(rworldmap::getMap("high"))
countriesArctic <- sf::st_intersection(countries,arctic)
#ice sheet list (vector names correspond to age)







#Load data
project <- 'RAW'
dir <- '/Users/chrishancock/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/Research/Manuscript/'
RAWts     <- readRDS(file.path(dir,project,'Data','RAW_LiPDts.RDS'))
shiftsAll <- readRDS(file.path(dir,project,'Data',paste(project,'meanShifts.RDS',sep='_')))

#Filter data
shifts <- shiftsAll %>% filter(interpVar %in% c('Temperature')) %>% filter(interpSeason %in% c('Annual','Summer','Winter'))
fTS <- RAWts[(pullTsVariable(RAWts,'paleoData_TSid') %in% shifts$paleoData_TSid)]
binvec <- sort(unique(c(shifts$time_start,shifts$time_end)))
binyrs <- sort(unique(c(shifts$time_mid)))

#find weights
dggs <- dgconstruct(res=4,pole_lat_deg=90,topology='HEXAGON')
polygons <- dggridR::dgearthgrid(dggs)
#locate each record's cell
shifts <- shifts %>% dplyr::mutate(cell=(dggridR::dgGEO_to_SEQNUM(dggs,geo_longitude,geo_latitude)$seqnum))
#add weights which consider number of record within each cell and the number of records at each site
shifts$weight <- NA
for (c in unique(shifts$cell)){
  cellShifts <- (shifts %>% filter (cell == c))
  cellweight <- 1/length(unique(cellShifts$dataSetName))
  for (s in unique(cellShifts$site)){
    proxyweight <- cellweight/length(unique((cellShifts %>% filter(site == s))$paleoData_TSid))
    shifts$weight[which((shifts$cell==c) & (shifts$site==s))] <- proxyweight
  }
}

# Add weight and climateInterpretation1_interpDirection (which is missing (or labeled differently for some records))
for (i in 1:length(fTS)){
  sel <- (shifts%>%filter(paleoData_TSid==fTS[[i]]$paleoData_TSid))
  fTS[[i]]$climateInterpretation1_interpDirection <- sel$interpDir[1]
  fTS[[i]]$geo_weight <- sel$weight[1]
}


#do composite
ensOut <- compositeEnsembles2(
  fTS = fTS,
  nens = 100,
  binvec = binvec,
  binFun = sampleEnsembleThenBinTs,
  ageVar = "ageEnsemble",
  uncVar = "paleoData_temperature12kUncertainty",
  weights = "geo_weight",
  stanFun = standardizeOverRandomInterval,
  searchRange = c(1000,11000),
  duration = 4000,
  normalizeVariance = TRUE,
  gaussianizeInput = TRUE,
  scale = FALSE,
  spread = TRUE,
  samplePct = 1,
)

plot(ensOut)
print(ensOut)













tsids <- unique(shifts$paleoData_TSid)

df <- data.frame()
pb <- progress_bar$new(total = length(tsids), format = "[:bar] :percent | ETA: :eta")
for (i in 1:length(tsids)){
  sel <- (shifts %>% filter(paleoData_TSid==tsids[i]))[1,28:length(names(shifts))]
  binmat <- matrix(NA,(length(binyrs)),ncol(sel$time[[1]]))
  for (c in 1:ncol(binmat)){
    binmat[,c] <- geoChronR::bin(time = sel$time[[1]][,c],values = sel$paleoData_values[[1]][,c],bin.vec = binvec)[,2]
  }
  if (sel$interpDir < 0){
    ages <- binyrs[unlist(apply(binmat,2,which.min))]
  }else{
    ages <- binyrs[unlist(apply(binmat,2,which.max))]
  }
  mode <- function(x){
    ux <- unique(x)
    m <- ux[which.max(tabulate(match(x, ux)))]
    return(m)
  }
  sel$agesMax <- list(ages)
  sel$agesMaxMode <- mode(ages)
  df <- rbind(df,sel)
  pb$tick()

}

plotMaxValAge <- function(df_raw,age.range,poly=NA, color.breaks = c(0.01,0.05,0.1,0.2,0.8,0.9,0.95,0.99),color.pal='RdBu'){
  #color.breaks <- c(0.001,0.01,0.025,0.05,0.95,0.975,0.99,0.999)
  #colorvec <- RColorBrewer::brewer.pal(length(age.range)+2, 'oranges')[2:(length(age.range)+2)]
  #colorvec[which(colorvec=="#F7F7F7")]<-'grey'
  # Point data
  points <- df_raw %>%  st_as_sf(coords = c("geo_longitude", "geo_latitude"), crs = 4326) #%>% arrange(desc(pvalue_either))
  # Convert pval and bin according to colorvec
  #points <- points[unlist(lapply(points$null_probability,sum)) > 0,]
  #points$pvalNet <- purrr::map2_dbl(points$pvalue_positive, points$pvalue_negative, handleNet)
  #points <- pval2factor(points,"pvalNet","pvalBinned",color.breaks)
  #
  cells <- data.frame()
  plotcells <- TRUE
  if (is.na(poly)){
    plotcells<-FALSE
  } else if (poly == 'cell'){
    cells <- data.frame(cell = sort(unique(points$cell)),ageSig=NA)
    for (i in 1:nrow(cells)){
      df_cell <- points %>% filter(cell==cells$cell[i])
      if (!("weight" %in% names(df_cell))){df_cell$weight<-NA}
      ages  <- sort(unique(df_cell$time_mid))
      pvals <- c() #
      counts <- c()
      for (age in ages){
        df_cell_age <- df_cell %>% filter(time_mid==age)
        out <- calculateMultiTestSignificance(df_cell_age,weights=df_cell_age$weight,n.ens=median(df_raw$null.hypothesis.n))
        pvals<- c(pvals, out$pvalEither) #handleNet(out$pvalPos,out$pvalNeg))
        counts <- c(counts, out$allEventPos)
      }
      #If a pvalie tie, go with the one with more positive shifts detected regardless of null
      if (min(pvals)<0.2){
        cells$ageSig[i] <- ages[pvals==min(pvals)][which.max(counts[pvals==min(pvals)])]
      }
      #cells$pvalNet[i] <- out$pvalNet
    }
    #
    cells <- merge(polygonsArctic,cells,by.x="seqnum",by.y="cell")
  } else(stop("poly must be 'cell' or 'region'"))
  #
  # Plot
  map <- ggplot() #+ actR_ggtheme()
  #Plot polygons
  if (plotcells){
    #cells$ageSig2 <- as.factor(cells$ageSig)

    #cells <- pval2factor(cells,"pvalNet","pvalBinned",color.breaks)
    map<-map+geom_sf(data=cells, aes(fill=ageSig2), color=NA)
  }
  map <- map + geom_sf(data=countriesArctic, fill=NA, color='grey30', linewidth=0.07)
  if (plotcells){
    map <- map + geom_sf(data=regions, fill=NA, color='black', linewidth=0.7)
  }
  map<-map +
    geom_sf(data=arcticCircle, fill=NA, color='black', linewidth=0.7) +
    geom_sf(data=points, aes(fill=as.factor(agesMaxMode)), color='black',shape=21,size=3) +
    #scale_fill_manual(values=colorvec, drop = FALSE,name='net\np-value',na.value = "white")+
    scale_fill_viridis_d(direction=-1)+
    coord_sf(crs =  sp::CRS("ESRI:102016"),ylim=c(-3300000,3300000),xlim=c(-3300000,3300000)) +
    #labs(title=paste0(min(age.range/1000),'-',max(age.range/1000),' ka \n age of most significant shift: '))+
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),        ## <- this line
      panel.border = element_rect(color=NA,fill=NA),
      panel.grid.major = element_blank(), #(color='black',linewidth=0.1),
      # panel.ontop = F
    )
  return(map)
}
plotMaxValAge(df)

z <- (df %>% filter(agesMaxMode<1000) )
i <- 3
hist(z[i,]$agesMax[[1]])
plot(apply(z[i,]$time[[1]],1,median),apply(z[i,]$paleoData_values[[1]],1,median))

