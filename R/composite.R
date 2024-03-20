#
#' Composite with uncertainty
#'
#' @param fTS TS object of the timeseries that you want to composite
#' @param binvec A vector of bin edges for binning
#' @param spread An option to spread values over adjacent empty years (if TRUE, spreadPaleoData() is used by binFun)
#' @param stanFun Function to use for standardization
#' @param uncVar Name of the uncertainty variable name in fTS. This is used for stanFun = sampleEnsembleThenBinTs but not simpleBinTs. If left NA or NULL, the default (1.5) will be applied.
#' @param ageVar Name of the time vector in fTS
#' @param gaussianizeInput Optionally gaussianize input before compositing
#' @param alignInterpDirection Optionally invert timeseries with negative interpretation_direction
#' @param scope interpretation_scope for alignInterpDirection (default = 'climate')
#' @param binFun Function to use for binning
#' @param ... parameters to pass to stanFun
#'
#'
#' @return a list. The main composite is in the list$composite matrix (needs a better definition here) TODO
#' @export
compositeEnsembles <- function(fTS,
                               binvec,
                               spread = TRUE,
                               stanFun = standardizeMeanIteratively,
                               ageVar = "age",
                               uncVar = NA,
                               gaussianizeInput = FALSE,
                               alignInterpDirection = TRUE,
                               scope = "climate",
                               binFun = sampleEnsembleThenBinTs,
                               ...){
  binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))
  #Align paleodata ages to common binstep
  if (is.character(uncVar)){
    binMatR <- as.matrix(purrr::map_dfc(fTS,binFun,binvec,ageVar = ageVar, uncVar=uncVar, spread = spread,gaussianizeInput = gaussianizeInput,alignInterpDirection = alignInterpDirection, scope = scope))
  }
  else{
    binMatR <- as.matrix(purrr::map_dfc(fTS,binFun,binvec,ageVar = ageVar,                spread = spread,gaussianizeInput = gaussianizeInput,alignInterpDirection = alignInterpDirection, scope = scope))
  }
  binMatR[is.nan(binMatR)] <- NA
  #Compute composites
  compMat <- stanFun(ages = binAges,pdm = binMatR,...)
  #Format output
  #which records contributed?
  gm <- is.finite(compMat)
  comp <- rowMeans(compMat,na.rm = TRUE)
  count <- colSums(gm)
  contributed <- which(count > 0)
  timeCount <- rowSums(gm)
  #Return output
  return(list(composite = comp, count = count, contributed = contributed,timeCount = timeCount))
}


#
#' Composite with uncertainty
#'
#' @param fTS TS object of the timeseries that you want to composite
#' @param binvec A vector of bin edges for binning
#' @param nens = The number of iterations to randomize over
#' @param spread An option to spread values over adjacent empty years (if TRUE, spreadPaleoData() is used by binFun)
#' @param stanFun Function to use for standardization
#' @param uncVar Uncertainty value. Enter a string to identify the variable name in fTS. Use a number to provide a unifomr value to all records. This is used for stanFun = sampleEnsembleThenBinTs but not simpleBinTs.
#' @param ageVar Name of the time vector in fTS
#' @param gaussianizeInput Optionally gaussianize input before compositing
#' @param alignInterpDirection Optionally invert timeseries with negative interpretation_direction
#' @param scope interpretation_scope for alignInterpDirection (default = 'climate')
#' @param binFun Function to use for binning
#' @param ... parameters to pass to stanFun
#'
#'
#' @return a list. The main composite is in the list$composite matrix (needs a better definition here) TODO
#' @export
compositeEnsembles2 <- function(fTS,
                                binvec,
                                nens = 100,
                                spread = TRUE,
                                stanFun = standardizeMeanIteratively,
                                ageVar = "age",
                                uncVar = NA,
                                gaussianizeInput = FALSE,
                                alignInterpDirection = TRUE,
                                scope = "climate",
                                binFun = sampleEnsembleThenBinTs,
                                ...){
  # Bin ages
  binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))
  # Do math
  doParallel::registerDoParallel(2)
  stanMatList <- foreach(i = 1:nens) %dopar% {
    set.seed(i)
    #Align paleodata ages to common binstep
    if (is.character(uncVar)){ #use uncVar in lipd file
      binMatR <- as.matrix(purrr::map_dfc(fTS,binFun,binvec,ageVar = ageVar, uncVar=uncVar,     spread = spread, gaussianizeInput = gaussianizeInput, alignInterpDirection = alignInterpDirection, scope = scope))
    } else if (is.numeric(uncVar)){ #use a default number
      binMatR <- as.matrix(purrr::map_dfc(fTS,binFun,binvec,ageVar = ageVar, defaultUnc=uncVar, spread = spread, gaussianizeInput = gaussianizeInput, alignInterpDirection = alignInterpDirection, scope = scope))
    } else{
      binMatR <- as.matrix(purrr::map_dfc(fTS,binFun,binvec,ageVar = ageVar,                    spread = spread, gaussianizeInput = gaussianizeInput, alignInterpDirection = alignInterpDirection, scope = scope))
    }
    #Standardize
    stanMatR <- stanFun(ages=binAges,pdm=binMatR,...)
    return(cbind(age=binAges,iteration=i,stanMatR))
  }
  # Summarize
  stan_df <- as.data.frame(abind::abind(stanMatList,along=1))
  # Composite means (each column is a composite mean of proxies)
  compMat <-stan_df %>%
    dplyr::transmute(age, iteration, mean = rowMeans(select(.,-c(age,iteration)), na.rm=T)) %>%
    dplyr::group_by(age,iteration) %>%
    tidyr::pivot_wider(names_from = iteration, values_from = mean) %>%
    dplyr::rename_at(vars(-age),~ paste0("comp", .))
  # Standardized proxy time series (each column is a proxy mean of iterations)
  proxyMat <- stan_df %>%
    select(-iteration) %>%
    group_by(age) %>%
    summarize(across(everything(), mean, na.rm=T)) %>%
    setNames(c('age',pullTsVariable(fTS,'paleoData_TSid')))
  # Count of proxy values used at each time
  proxyPctUsed <-data.frame(age=stan_df$age,!is.na(stan_df%>%select(.,-c(age,iteration)))) %>%
    group_by(age) %>%
    summarize(across(everything(), sum)/nens) %>%
    setNames(c('age',pullTsVariable(fTS,'paleoData_TSid')))
  # Return
  return(list(composite = compMat, proxyVals = proxyMat, proxyUsed = proxyPctUsed))
}

#
#' Plot compositeEnsembles output
#'
#' @param compositeList The output of compositeEnsembles()
#' @param ageUnits units for the age variable
#' @param valUnits units for the composite
#' @param title a title for the combined plot (not used if combine == FALSE)
#' @param combine output a combined figure using egg::ggarrange or provide a list of subplots
#' @param ... parameters to pass to geoChronR::plotTimeseriesEnsRibbons
#'
#' @return a combined or list of gg plots
#' @export
plotComposite <- function(compositeList, ageUnits = 'yr BP',valUnits='standardized',title='composite proxy record',combine=TRUE,...){
  #Plot
  ensRibbon <- geoChronR::plotTimeseriesEnsRibbons(X=list(values = compositeList$composite$age, units=ageUnits, variableName='age'),
                                        Y=list(values = compositeList$composite[,-1],units=valUnits, variableName='composite anomaly'),
                                        ...)+theme_bw()+
                                        scale_x_reverse(limits=rev(ar),expand=c(0,0),name=paste0('age (',ageUnits,')'))
  if (combine){ensRibbon <- ensRibbon + scale_x_reverse(limits=rev(ar),expand=c(0,0),name=' ',position="top")}
  #
  #Plot Count
  bs = abs(median(diff(compositeList$composite$age)))
  ar = c(min(compositeList$proxyVals,na.rm=T)-bs/2,max(compositeList$proxyVals,na.rm=T)+bs/2)
  tsAvailability <- ggplot2::ggplot(compositeList$proxyUsed)+
    geom_bar(stat="identity",aes(x=age,y=apply(compositeList$proxyUsed %>% select(-age),1,sum)),width=bs)+
    scale_y_continuous(limits=c(0,ncol(compositeList$proxyUsed)),expand=c(0,0),name='count')+
    scale_x_reverse(limits=rev(ar),expand=c(0,0),name=paste0('age (',ageUnits,')')) + theme_bw()
  #
  #Combine
  plots <- list(compositeEns=ensRibbon,dataUsed=tsAvailability)
  if (combine){
    plots <-   egg::ggarrange(plots=plots,ncol=1,heights=c(2,1),top=title)
  }
  return(plots)
}





#TO DO
#Add spatial composites ?
#How to sample paleo uncertainty? sampleTHenBin

#'
#' compositeSCC (Standard Calibrated Composite) with presets for Temp12k
#'
#' @param fTS TS object of the timeseries that you want to composite
#' @param binvec A vector of bin edges for binning
#' @param interval The ager range used to create the standardization interval
#' @param ... parameters to pass to stanFun
#'
#' @return a list. The main composite is in the list$composite matrix (needs a better definition here) TODO
#' @export
#'
compositeSCC <- function(fTS, binvec, interval=c(3000,5000), ...){
  results <- compositeEnsembles(fTS,
                                binvec,
                                stanFun = standardizeOverInterval,
                                interval = interval,
                                spread = TRUE,
                                gaussianizeInput = FALSE,
                                normalizeVariance = FALSE,
                                ...
                                )
  return(results)
}


#'
#' compositeDCC (Dynamic Calibrated Composite) with presets for Temp12k
#'
#' @param fTS TS object of the timeseries that you want to composite
#' @param binvec A vector of bin edges for binning
#' @param duration The length of time used to create the standardization interval
#' @param searchRange An age range from which duration is sampled within
#' @param spread An option to spread values over adjacent empty years (if TRUE, spreadPaleoData() is used by binFun)
#' @param gaussianizeInput Optionally gaussianize input before compositing
#' @param alignInterpDirection Optionally invert timeseries with negative interpretation_direction
#' @param ... parameters to pass to stanFun
#'
#' @return a list. The main composite is in the list$composite matrix (needs a better definition here) TODO
#' @export
#'
compositeDCC <- function(fTS, binvec, duration=3000, searchRange=c(0,7000),spread=TRUE,gaussianizeInput=FALSE, alignInterpDirection=TRUE,...){
  results <- compositeEnsembles(fTS,
                                binvec,
                                stanFun     = standardizeMeanIteratively,
                                duration    = duration,
                                searchRange = searchRange,
                                spread            = spread,
                                gaussianizeInput  = gaussianizeInput,
                                alignInterpDirection = alignInterpDirection,
                                ...
                                )
  return(results)
}

compositeGAM <- function(GAM){
  print('This method is closely related to SCC, but instead of computing the mean of the records within every 100-year interval, it fits a generalized additive model (GAM) through the ensemble. ')
}

compositeCPS <- function(CPS){
  print('This method is closely related to SCC, but instead of computing the mean of the records within every 100-year interval, it fits a generalized additive model (GAM) through the ensemble. ')
}

compositePAI <- function(PAI){
  print('This method is closely related to SCC, but instead of computing the mean of the records within every 100-year interval, it fits a generalized additive model (GAM) through the ensemble. ')
}
# library(sf)
# library(dplyr)
# library(ggplot2)
# binLatLonGrid <- function(lon_range=c(-180,180),lat_range=c(-90,90),cellsize=c(360,30)){
#   lon_vec <- seq(lon_range[1], lon_range[2], cellsize[1])
#   lat_vec <- seq(lat_range[1], lat_range[2], tail(cellsize, n=1))
#   #
#   coords <- as.data.frame(expand.grid(lon = lon_vec, lat = lat_vec))
#   colnames(coords) <- c('lon','lat')
#   #
#   polygons <- as.data.frame(coords) %>%
#     sf::st_as_sf(coords=c('lon','lat'), crs=4326) %>%
#     sf::st_make_grid(what="polygons", cellsize=cellsize) %>%
#     sf::st_as_sf() %>%
#     mutate(lon=sf::st_coordinates(sf::st_centroid(.))[,1]) %>%
#     mutate(lat=sf::st_coordinates(sf::st_centroid(.))[,2]) %>%
#     mutate(area_km2 = units::drop_units(sf::st_area(.)/1000000))
#   return(polygons)
# }
#
# countries <- st_as_sf(rworldmap::getMap("high"))
#
#
# ggplot(binLatLonGrid(cellsize=c(360,30)))+
#   geom_sf(color='red')+
#   geom_sf(data=countries,fill=NA)+
#   coord_sf()
#
# ggplot(binLatLonGrid(cellsize=c(10,5)))+
#   geom_sf(color='red')+
#   geom_sf(data=countries,fill=NA)+
#   coord_sf()
#
#
# binEqualAreaGrid <- function(res=3,metric=FALSE, resround='down'){
#   dggs <- dggridR::dgconstruct(res=res,topology='HEXAGON',pole_lat_deg = 90,pole_lon_deg = 180)
#   polygons <- dggridR::dgearthgrid(dggs) %>%
#     dplyr::mutate(lon=sf::st_coordinates(sf::st_centroid(.))[,1]) %>%
#     dplyr::mutate(lat=sf::st_coordinates(sf::st_centroid(.))[,2]) %>%
#     dplyr::mutate(area_km2 = units::drop_units(sf::st_area(.)/1000000))
#   return(polygons)
# }
#
# z<-binEqualAreaGrid(res=4)
# z<-st_wrap_dateline(z, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)
# library(dggridR)
#
# #Construct a global grid with cells approximately 1000 miles across
# dggs          <- dgconstruct(spacing=1000, metric=FALSE, resround='down')
#
# #Load included test data set
# data(dgquakes)
#
# #Get the corresponding grid cells for each earthquake epicenter (lat-long pair)
# dgquakes$cell <- dgGEO_to_SEQNUM(dggs,dgquakes$lon,dgquakes$lat)$seqnum
#
# #Converting SEQNUM to GEO gives the center coordinates of the cells
# cellcenters   <- dgSEQNUM_to_GEO(dggs,dgquakes$cell)
#
# #Get the number of earthquakes in each cell
# quakecounts   <- dgquakes %>% group_by(cell) %>% summarise(count=n())
#
# #Get the grid cell boundaries for cells which had quakes
# grid          <- dgcellstogrid(dggs,quakecounts$cell)
#
# #Update the grid cells' properties to include the number of earthquakes
# #in each cell
# grid          <- merge(grid,quakecounts,by.x="seqnum",by.y="cell")
#
# #Make adjustments so the output is more visually interesting
# grid$count    <- log(grid$count)
# cutoff        <- quantile(grid$count,0.9)
# grid          <- grid %>% mutate(count=ifelse(count>cutoff,cutoff,count))
#
# #Get polygons for each country of the world
# countries <- map_data("world")
# sf_use_s2(FALSE)
# wrapped_grid = st_wrap_dateline(grid, options = c("WRAPDATELINE=YES","DATELINEOFFSET=0"))
#
# ggplot() +
#   geom_sf(data=(grid %>% filter(sf::st_is_valid(.))), aes(fill=count), color=alpha("white", 0.4)) +
#   geom_sf(data=countries, fill=NA, color="black")   +
#   scale_fill_gradient(low="blue", high="red")#+coord_sf(crs="+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83")
#
#
#
# #Get polygons for each country of the world
# countries <- map_data("world")
# ggplot()+
#   geom_sf(data =z,color='red')+
#   geom_sf(data=countries,fill=NA)+
#   coord_sf(crs="+proj=laea")
# #+
#   #
#   c#oord_sf(crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +datum=WGS84 +units=m +no_defs")
#
#
#
#
#
#
#
#
#
#
# ggplot(binIPCCRegions(4))+
#   geom_sf(color='red')+
#   geom_sf(data=countries,fill=NA)+
#   coord_sf()
#
#
# binIPCCRegions <- function(x){
#   refregions <- load(url('https://github.com/SantanderMetGroup/ATLAS/blob/main/reference-regions/IPCC-WGI-reference-regions-v4_R.rda?raw=true'), verbose = F)
#   polygons <- sf::st_as_sf(as(IPCC_WGI_reference_regions_v4, "SpatialPolygons"))   %>%
#     mutate(lon=sf::st_coordinates(sf::st_centroid(.))[,1]) %>%
#     mutate(lat=sf::st_coordinates(sf::st_centroid(.))[,2]) %>%
#     mutate(area_km2 = units::drop_units(sf::st_area(.)/1000000)) %>%
#     mutate(region = names(IPCC_WGI_reference_regions_v4@polygons))
#   return(polygons)
# }
