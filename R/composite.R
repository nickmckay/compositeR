#' Convert timeseries input to a lipdR TS list
#'
#' compositeR functions work on a list with one element per timeseries. This
#' accepts that, or a lipdR timeseries tibble (`lipdR::ts2tibble()`), or any
#' data frame with one row per timeseries.
#'
#' @param fTS a TS list, lipd_ts_tibble, or data frame
#' @return a TS list
#' @keywords internal
asTsList <- function(fTS){
  if(is.data.frame(fTS) || inherits(fTS, "lipd_ts")){
    return(lipdR::as.lipdTs(fTS))
  }
  fTS
}

#'
#' Composite ensembles
#'
#' Create paleoclimate composite ensembles from a set of LiPD timeseries. Written by Chris Hancock (2024) as compositeEnsembles2(); since compositeR 0.4.0 it is the main compositing function.
#'
#' @import dplyr
#'
#' @param fTS the timeseries to composite: a lipdR TS list (from `lipdR::extractTs()`), or a timeseries tibble/data frame with one row per timeseries (from `lipdR::ts2tibble()`, or any data frame with list columns for the age and paleoData_values)
#' @param binvec A vector that describes the edges of the bins to bin to
#' @param nens the number of ensembles to create
#' @param stanFun function to use for standardization (either standardizeOverRandomInterval, standardizeMeanIteratively, or standardizeOverInterval)
#' @param binFun function to use for binning (either sampleEnsembleThenBinTs or simpleBinTS)
#' @param uncVar uncertainty variable  \cr - Use a string to identify a fTS variable name or use a number to provide a uniform value to all records.  \cr - This is used if binFun == sampleEnsembleThenBinTs but not for simpleBinTs.
#' @param uncUnits units of a numeric `uncVar`: "native" (default) for each record's own units, or "sd" for multiples of each record's standard deviation (e.g. `uncVar = 0.5, uncUnits = "sd"` gives every record noise with half its SD). Must be "native" when `uncVar` names a variable, since those are measured uncertainties.
#' @param weights weights for calculating the composite mean. \cr - For example, you may want to weight each record in a region relative to the number of total records from a single site.  \cr - Use a string to identify a fTS variable name, otherwise all values will be weighted equally.
#' @param samplePct the percent of records to used for each ensemble.  \cr - If length(fTS)*samplePct <= 3, this is ignored.  \cr - Default = 1 which means that 100 percent of records are used for each ensemble
#' @param scale scale each composite ensemble member to have a mean of zero and unit variance at the end of the compiting
#' @param ageVar the time variable to composite on. The default, "time", uses each record's time axis converted to `timeUnits` with [lipdR::getTsTime()], so records in year CE, yr BP or ka can be composited together, and an age ensemble (`ageEnsemble`) is used when a record has one. Any other name (e.g. "age" or "ageEnsemble") is used as it is, without conversion, as before compositeR 0.5.0.
#' @param timeUnits the time scale to composite on when `ageVar = "time"` (default "yr BP"). `binvec` is in these units. Also used to label the output.
#' @param spread "Spread" data before binning to prevent aliasing of bins (default = TRUE)
#' @param gaussianizeInput "Gaussianize" the data before compositing to avoid impacts of irregularly distributed data (default = TRUE)
#' @param alignInterpDirection multiply each record by -1 if its interpretation direction is negative, so that all records point the same way (see [interpDirectionSign()])
#' @param scope which interpretation scope to align with (default "climate"). Interpretations without a scope are used if none match. NA uses the first interpretation with a direction.
#' @param ... additional options to pass to stanFun
#' @param verbose should a progress bar be printed?
#'
#' @return a composite list class with three tibbles:
#' \cr (1) The composite [dims = length(binvec) x (nens + 1 for age column)]
#' \cr (2) the standardized values [dims = length(binvec) x (length(fTS) + 1 for age column)], and
#' \cr (3) a count of when proxy values are used in the compilation [same dims as 2]
#' @export
compositeEnsembles <- function(fTS,
                               binvec,
                               nens = 100,
                               stanFun = standardizeOverRandomInterval,
                               ageVar = "time",
                               timeUnits = "yr BP",
                               uncVar = NA,
                               uncUnits = c("native", "sd"),
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
  fTS <- asTsList(fTS)
  uncUnits <- match.arg(uncUnits)
  if(is.character(uncVar) && uncUnits == "sd"){
    stop("uncUnits = 'sd' applies to a numeric uncVar. Uncertainties read from the timeseries ('", uncVar, "') are in the records' own units.")
  }

  # Put every record on the same time scale # #############################
  binAgeVar <- ageVar
  if(identical(ageVar, "time")){
    fTS <- addCompositeTime(fTS, timeUnits, ensemble = !identical(binFun, simpleBinTs))
    binAgeVar <- ".compositeTime"
  }

  # Bin ages # #############################
  if(!is.numeric(binvec)){stop("binvec must be numeric")}
  binvec <- sort(binvec)
  binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

  # Warnings # #############################
  # Check that necessary variables are provided in lipd files
  for (name in c("paleoData_TSid",binAgeVar)){
    check <- unlist(lapply(lipdR::pullTsVariable(fTS,name),is.null))
    if (sum(check)>0){ stop(paste0("ERROR: No '",name,"' variable for fTS[c(",(paste(which(!check),collapse=', ')),")]")) }
  }
  # binvec
  if(sum(stats::median(diff(binvec)) != diff(binvec)) > 0){warning("binvec has uneven spacing. This vector should be created using seq(ageMin,ageMax,binsize")}
  # alignInterpDirection
  if (alignInterpDirection){
    signs <- vapply(fTS, interpDirectionSign, numeric(1), scope = scope)
    if(any(is.na(signs))){
      warning(paste0("alignInterpDirection: no usable ", ifelse(is.na(scope), "", paste0("'", scope, "' ")),
                     "interpretation direction for fTS[c(", paste(which(is.na(signs)), collapse = ", "),
                     ")]. These records are not flipped."))
    }
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
        fTS, binFun, binvec, ageVar=binAgeVar, uncVar=uncVar, #uncVar=uncVar
        spread=spread, gaussianizeInput=gaussianizeInput, alignInterpDirection=alignInterpDirection, scope=scope
        )))
    } else if (is.numeric(uncVar)){ #use a default number
      binMatR <- suppressMessages(as.matrix(purrr::map_dfc(
        fTS, binFun, binvec, ageVar=binAgeVar, defaultUnc=uncVar, uncUnits=uncUnits,
        spread=spread, gaussianizeInput=gaussianizeInput, alignInterpDirection=alignInterpDirection, scope=scope
        )))
    } else{
      binMatR <- suppressMessages(as.matrix(purrr::map_dfc(
        fTS, binFun, binvec, ageVar=binAgeVar,
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
    dplyr::summarize(dplyr::across(dplyr::everything(), \(x) mean(x, na.rm = TRUE))) %>%
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
              proxyUsed = proxyPctUsed[,-1],
              timeUnits = timeUnits)
  out[['parameters']] <- data.frame(nens=nens,
                                    ageVar=ageVar,
                                    timeUnits=timeUnits,
                                    uncVar=uncVar,
                                    weights=weights,
                                    samplePct=samplePct,
                                    scale=scale,
                                    spread=spread,
                                    gaussianizeInput=gaussianizeInput,
                                    alignInterpDirection=alignInterpDirection,
                                    scope=scope)
  #`paramaters` (misspelled) is retained for one deprecation cycle; use
  #`parameters`.
  out[['paramaters']] <- out[['parameters']]
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


#' compositeEnsembles2 (deprecated)
#'
#' Renamed to [compositeEnsembles()] in compositeR 0.4.0. This alias calls it
#' with the same arguments, and will be removed in a future release.
#'
#' @param ... arguments passed to [compositeEnsembles()]
#' @return see [compositeEnsembles()]
#' @export
compositeEnsembles2 <- function(...){
  if(is.null(getOption("compositeR.compositeEnsembles2.deprecated"))){
    warning("compositeEnsembles2() has been renamed compositeEnsembles(). The old name will be removed in a future release.")
    options(compositeR.compositeEnsembles2.deprecated = TRUE)
  }
  compositeEnsembles(...)
}


#' Add each record's time, on one scale, for compositing
#'
#' @param fTS a TS list
#' @param timeUnits the target time scale
#' @param ensemble use age ensembles where records have them?
#' @return fTS with a `.compositeTime` field on every record
#' @keywords internal
addCompositeTime <- function(fTS, timeUnits, ensemble = TRUE){
  for(i in seq_along(fTS)){
    tm <- tryCatch(lipdR::getTsTime(fTS[[i]], units = timeUnits, ensemble = ensemble),
                   error = function(e){
                     id <- fTS[[i]]$paleoData_TSid %||% fTS[[i]]$dataSetName %||% paste0("fTS[[", i, "]]")
                     stop("Can't put ", id, " on the ", timeUnits, " scale: ", conditionMessage(e), call. = FALSE)
                   })
    attr(tm, "timeUnits") <- NULL
    fTS[[i]]$.compositeTime <- tm
  }
  fTS
}

`%||%` <- function(x, y){
  if(is.null(x)) y else x
}

# Does a time scale count forward (year CE) rather than back (yr BP, ka)?
isProgradeUnits <- function(units){
  !is.null(units) && !is.na(units) &&
    grepl("\\b(AD|CE)\\b", units, ignore.case = TRUE) &&
    !grepl("\\bBC", units, ignore.case = TRUE)
}
