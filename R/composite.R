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
  binMatR <- as.matrix(purrr::map_dfc(fTS,
                                      binFun,
                                      binvec,
                                      ageVar = ageVar,
                                      spread = spread,
                                      gaussianizeInput = gaussianizeInput,
                                      alignInterpDirection = alignInterpDirection,
                                      scope = scope))
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
#'
#' @param fTS TS list of the timeseries that you want to composite
#' @param binvec A vector that describes the edges of the bins to bin to
#' @param nens the number of ensembles to create
#' @param stanFun function to use for standardization (either standardizeOverRandomInterval, standardizeMeanIteratively, or standardizeOverInterval)
#' @param binFun function to use for binning (either sampleEnsembleThenBinTs or simpleBinTS)
#' @param uncVar uncertainty variable  \cr - Use a string to identify a fTS variable name or use a number to provide a uniform value to all records.  \cr - This is used if stanFun == sampleEnsembleThenBinTs but not for simpleBinTs.
#' @param weights weights for calculating the composite mean. \cr - For example, you may want to weight each record in a region relative to the number of total records from a single site.  \cr - Use a string to identify a fTS variable name, otherwise all values will be weighted equally.
#' @param samplePct the percent of records to used for each ensemble.  \cr - If length(fTS)*samplePct <= 3, this is ignored.  \cr - Default = 1 which means that 100 percent of records are used for each ensemble
#' @param scale scale each composite ensemble member to have a mean of zero and unit variance at the end of the compiting
#' @param ageVar age variable to use
#' @param spread "Spread" data before binning to prevent aliasing of bins (default = TRUE)
#' @param gaussianizeInput "Gaussianize" the data before compositing to avoid impacts of irregularly distributed data (default = TRUE)
#' @param alignInterpDirection align data by interpretation direction
#' @param scope which scope of interpretation to use
#' @param ... additional options to pass to stanFun
#' @param verbose should a progress bar be printed?
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
                                verbose = TRUE,
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
  print(glue::glue("binning and aligning the data. This may take some time depending on the number of records (n = {length(fTS)}) and number of ensembles (n = {nens})"))
  if(verbose){pb <- progress::progress_bar$new(total = nens, format = "[:bar] :percent | ETA: :eta")}
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
    sampleMin <- min(3,length(fTS)) #If number of records * samplePct <= sampleMin, samplePct will be ignored
    binMatR[,-(sample(ncol(binMatR),max(ceiling(ncol(binMatR)*samplePct),sampleMin)))] <- NA #TODO fix warnings from this
    #Standardize
    stanMat <- stanFun(ages=binAges, pdm=binMatR, ...)
    #Return
    stanMatList[[i]] <- cbind(age=binAges, iteration=i, stanMat)
    if(verbose){pb$tick()}
  }


  #TODO: standardize(document)

  # Summarize
  stan_df <- as.data.frame(abind::abind(stanMatList,along=1))
  # Composite means (each column is a composite mean of proxies)

  scaleNaRm <- function(x){
    out <- (x - mean(x,na.rm = TRUE)) / sd(x,na.rm = TRUE)
    return(out)
  }

  compMat <- stan_df %>%
    dplyr::transmute(age, iteration, mean = apply(dplyr::select(.,-c(age,iteration)),1,stats::weighted.mean,w = w, na.rm = TRUE)) %>%
    dplyr::group_by(age,iteration) %>%
    tidyr::pivot_wider(names_from = iteration, values_from = mean) %>%
    dplyr::rename_at(dplyr::vars(-age),~ paste0("comp", .))

  if(scale){
    compMat <- mutate(ungroup(compMat), across(starts_with("comp"),.fns = scaleNaRm))
  }

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
