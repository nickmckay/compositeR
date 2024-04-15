#
#' Composite with uncertainty
#'
#' @param fTS TS object of the timeseries that you want to composite
#' @param binvec A vector of bin edges for binning
#' @param spread An option to spread values over adjacent empty years (if TRUE, spreadPaleoData() is used by binFun)
#' @param stanFun Function to use for standardization
#' @param ageVar Name of the time vector in fTS
#' @param gaussianizeInput Optionally gaussianize input before compositing
#' @param alignInterpDirection Optionally invert timeseries with negative interpretation_direction
#' @param scope interpretation_scope for alignInterpDirection (default = 'climate')
#' @param binFun Function to use for binning
#' @param ... parameters to pass to stanFun
#'
#' @importFrom magrittr %>%
#'
#' @return a list. The main composite is in the list$composite matrix (needs a better definition here) TODO
#' @export
compositeEnsembles <- function(fTS,
                               binvec,
                               spread = TRUE,
                               stanFun = standardizeMeanIteratively,
                               ageVar = "age",
                               gaussianizeInput = FALSE,
                               alignInterpDirection = TRUE,
                               scope = "climate",
                               binFun = sampleEnsembleThenBinTs,
                               ...){

  binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

  #Align paleodata ages to common binstep
  binMatR <- as.matrix(purrr::map_dfc(fTS,binFun,binvec,ageVar = ageVar,spread = spread,gaussianizeInput = gaussianizeInput,alignInterpDirection = alignInterpDirection, scope = scope))
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



#'
#' create paleoclimate composite ensembles from a list of LiPD timeseries. A revised version of compositeEnsembles by Chris Hancock (2024)
#'
#' @import dplyr
#' @param fTS TS list of the timeseries that you want to composite
#' @param nens the number of ensembles to create
#' @param stanFun function to use for standardization (either standardizeOverRandomInterval, standardizeMeanIteratively, or standardizeOverInterval)
#' @param binFun function to use for binning (either sampleEnsembleThenBinTs or simpleBinTS)
#' @param uncVar uncertainty variable  \cr - Use a string to identify a fTS variable name or use a number to provide a uniform value to all records.  \cr - This is used if stanFun == sampleEnsembleThenBinTs but not for simpleBinTs.
#' @param weights weights for calculating the composite mean. \cr - For example, you may want to weight each record in a region relative to the number of total records from a single site.  \cr - Use a string to identify a fTS variable name, otherwise all values will be weighted equally.
#' @param samplePct the percent of records to used for each ensemble.  \cr - If length(fTS)*samplePct <= 3, this is ignored.  \cr - Default = 1 which means that 100 percent of records are used for each ensemble
#' @param scale TODO
#' @inheritParams sampleEnsembleThenBinTs
#' @inheritDotParams standardizeOverRandomInterval duration searchRange minN normalizeVariance
#'
#' @return a composite list class with three tibbles:
#' \cr (1) The composite [dims = length(binvec) x (nens + 1 for age column)]
#' \cr (2) the standardized values [dims = length(binvec) x (length(fTS) + 1 for age column)], and
#' \cr (3) a count of when proxy values are used in the compilation [same dims as 2]
#' @export
compositeEnsembles2 <- function(fTS,
                                binvec,
                                nens = 100,
                                stanFun = standardizeOverRandomInterval,
                                ageVar = "age",
                                uncVar = NA,
                                weights = NA,
                                samplePct = 1,
                                scale = FALSE,
                                spread = TRUE,
                                gaussianizeInput = FALSE,
                                alignInterpDirection = TRUE,
                                scope = "climate",
                                binFun = sampleEnsembleThenBinTs,
                                ...){
  # Bin ages # #############################
  if(!is.numeric(binvec)){stop("binvec must be numeric")}
  binvec <- sort(binvec)
  binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

  # Warnings # #############################
  # Check that necessary variables are provided in lipd files
  for (name in c("paleoData_TSid",ageVar)){
    check <- unlist(lapply(lipdR::pullTsVariable(fTS,name),is.null))
    if (sum(check)>0){ stop(paste0("ERROR: No '",name,"' variable for fTS[c(",(paste(which(!check),collapse=', ')),")]")) }
  }
  # binvec
  if(sum(stats::median(diff(binvec)) != diff(binvec)) > 0){warning("binvec has uneven spacing. This vector should be created using seq(ageMin,ageMax,binsize")}
  # alignInter
  if (alignInterpDirection){
    likelyname <- paste0(scope,"Interpretation1_interpDirection")
    if(sum(is.na(lipdR::pullTsVariable(fTS,likelyname)))>1){warning(paste0("Warning with alignInterpDirection and scope. Variable name '",likelyname,"' not found for fTS[ c(",paste(which(is.na(lipdR::pullTsVariable(fTS,likelyname))),collapse=', '),")]"))}
    if (sum(!(tolower(unique(lipdR::pullTsVariable(fTS,likelyname))) %in% c("positive","negative",-1,1,NA)))>0){ warning(paste0("Unexpected ",likelyname," found in fTS. This value should be positive/negative or 1/-1"))}
  }

  # weights
  if (is.na(weights)){ # weight all equally
    w <- rep(1,length(fTS))
  } else{
    w <- as.vector(lipdR::pullTsVariable(fTS,weights))
    if (is.numeric(w)){ # use this
      w[is.na(w)] <- max(w,na.rm=T) # if this value is missing, assign it with the max w value
    } else{ # if a character vector such as geo_SiteName, geo_Region, or archiveType, assign based on number of unique values
      for (i in unique(w)){
        w[w==i] <- 1/sum(w==i)
      }
      w <- as.numeric(w)
    }
  }

  # Bin and standardize the timeseries for each iteration  # #############################
  stanMatList <- list()
  print("binning and aligning the data. This may take some time depending on the number of records (length(fTS)) and number of ensembles (nens)")
  pb <- progress::progress_bar$new(total = nens, format = "[:bar] :percent | ETA: :eta")
  for (i in 1:nens){
    set.seed(i)
    #Align paleodata ages to common binstep
    if (is.character(uncVar)){ #use uncVar in lipd file
      binMatR <- suppressMessages(as.matrix(purrr::map_dfc(
        fTS, binFun, binvec, ageVar=ageVar, uncVar=uncVar, #uncVar=uncVar
        spread=spread, gaussianizeInput=gaussianizeInput, alignInterpDirection=alignInterpDirection, scope=scope
        )))
    } else if (is.numeric(uncVar)){ #use a default number
      binMatR <- suppressMessages(as.matrix(purrr::map_dfc(
        fTS, binFun, binvec, ageVar=ageVar, defaultUnc=uncVar, #defaultUnc=uncVar
        spread=spread, gaussianizeInput=gaussianizeInput, alignInterpDirection=alignInterpDirection, scope=scope
        )))
    } else{
      binMatR <- suppressMessages(as.matrix(purrr::map_dfc(
        fTS, binFun, binvec, ageVar=ageVar,
        spread=spread, gaussianizeInput=gaussianizeInput, alignInterpDirection=alignInterpDirection, scope=scope
        )))
    }
    #Sample different combination of records for iterations of the the composite
    sampleMin <- 3 #If number of records * samplePct <= sampleMin, samplePct will be ignored
    binMatR[,-(sample(ncol(binMatR),max(ceiling(ncol(binMatR)*samplePct),sampleMin)))] <- NA #TODO fix warnings from this
    #Standardize
    stanMat <- stanFun(ages=binAges, pdm=binMatR, ...)
    #Return
    stanMatList[[i]] <- cbind(age=binAges, iteration=i, stanMat)
    pb$tick()
  }

  #Scale
  if (scale){warning('scale not implemented yet. no scaling applied')}
  #TODO: standardize(document)

  # Summarize
  stan_df <- as.data.frame(abind::abind(stanMatList,along=1))
  # Composite means (each column is a composite mean of proxies)
  compMat <- stan_df %>%
    dplyr::transmute(age, iteration, mean = apply(dplyr::select(.,-c(age,iteration)),1,stats::weighted.mean, w=w, na.rm=T)) %>%
    dplyr::group_by(age,iteration) %>%
    tidyr::pivot_wider(names_from = iteration, values_from = mean) %>%
    dplyr::rename_at(dplyr::vars(-age),~ paste0("comp", .))
  # Standardized proxy time series (each column is a proxy mean of iterations)
  proxyMat <- stan_df %>%
    dplyr::select(-iteration) %>%
    dplyr::group_by(age) %>%
    dplyr::summarize(dplyr::across(dplyr::everything(), mean, na.rm=T)) %>%
    stats::setNames(c('age',lipdR::pullTsVariable(fTS,'paleoData_TSid')))
  # Count of proxy values used at each time
  proxyPctUsed <-data.frame(age=stan_df$age,!is.na(stan_df %>% dplyr::select(.,-c(age,iteration)))) %>%
    dplyr::group_by(age) %>%
    dplyr::summarize(dplyr::across(dplyr::everything(), sum)/nens) %>%
    stats::setNames(c('age',lipdR::pullTsVariable(fTS,'paleoData_TSid')))
  # Return
  out <- list(ages =  compMat$age,
              composite = compMat[,-1],
              proxyVals = proxyMat[,-1],
              proxyUsed = proxyPctUsed[,-1])
  out[['paramaters']] <- data.frame(nens=nens,
                                    ageVar=ageVar,
                                    uncVar=uncVar,
                                    weights=weights,
                                    samplePct=samplePct,
                                    scale=scale,
                                    spread=spread,
                                    gaussianizeInput=gaussianizeInput,
                                    alignInterpDirection=alignInterpDirection,
                                    scope=scope)
  return(new_composite(out))
}



#' Make S3 class for compositeEnsembles2
#'
#' @param x list output from compositeEnsembles2()
#'
#' @return list output from compositeEnsembles2() as S3 class "composite"
#' @export
new_composite <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x,class = c("paleoComposite",class(list())))
}



#'
#'
#' #TO DO
#' #Add spatial composites ?
#' #How to sample paleo uncertainty? sampleTHenBin
#'
#' #'
#' #' compositeSCC (Standard Calibrated Composite) with presets for Temp12k
#' #'
#' #' @param fTS TS object of the timeseries that you want to composite
#' #' @param binvec A vector of bin edges for binning
#' #' @param interval The ager range used to create the standardization interval
#' #' @param ... parameters to pass to stanFun
#' #'
#' #' @return a list. The main composite is in the list$composite matrix (needs a better definition here) TODO
#' #' @export
#' #'
#' compositeSCC <- function(fTS, binvec, interval=c(3000,5000), ...){
#'   results <- compositeEnsembles(fTS,
#'                                 binvec,
#'                                 stanFun = standardizeOverInterval,
#'                                 interval = interval,
#'                                 spread = TRUE,
#'                                 gaussianizeInput = FALSE,
#'                                 normalizeVariance = FALSE,
#'                                 ...
#'                                 )
#'   return(results)
#' }
#'
#'
#' #'
#' #' compositeDCC (Dynamic Calibrated Composite) with presets for Temp12k
#' #'
#' #' @param fTS TS object of the timeseries that you want to composite
#' #' @param binvec A vector of bin edges for binning
#' #' @param duration The length of time used to create the standardization interval
#' #' @param searchRange An age range from which duration is sampled within
#' #' @param spread An option to spread values over adjacent empty years (if TRUE, spreadPaleoData() is used by binFun)
#' #' @param gaussianizeInput Optionally gaussianize input before compositing
#' #' @param alignInterpDirection Optionally invert timeseries with negative interpretation_direction
#' #' @param ... parameters to pass to stanFun
#' #'
#' #' @return a list. The main composite is in the list$composite matrix (needs a better definition here) TODO
#' #' @export
#' #'
#' compositeDCC <- function(fTS, binvec, duration=3000, searchRange=c(0,7000),spread=TRUE,gaussianizeInput=FALSE, alignInterpDirection=TRUE,...){
#'   results <- compositeEnsembles(fTS,
#'                                 binvec,
#'                                 stanFun     = standardizeMeanIteratively,
#'                                 duration    = duration,
#'                                 searchRange = searchRange,
#'                                 spread            = spread,
#'                                 gaussianizeInput  = gaussianizeInput,
#'                                 alignInterpDirection = alignInterpDirection,
#'                                 ...
#'                                 )
#'   return(results)
#' }
#'
#' compositeGAM <- function(GAM){
#'   print('This method is closely related to SCC, but instead of computing the mean of the records within every 100-year interval, it fits a generalized additive model (GAM) through the ensemble. ')
#' }
#'
#' compositeCPS <- function(CPS){
#'   print('This method is closely related to SCC, but instead of computing the mean of the records within every 100-year interval, it fits a generalized additive model (GAM) through the ensemble. ')
#' }
#'
#' compositePAI <- function(PAI){
#'   print('This method is closely related to SCC, but instead of computing the mean of the records within every 100-year interval, it fits a generalized additive model (GAM) through the ensemble. ')
#' }
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
