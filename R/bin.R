#' removeConsecutiveDuplicates
#'
#' Remove consecutive duplicates within a time series.
#'
#' @param spreadValues a vector of binned paleoData_values
#' @return a matrix vector
#' @export
#'
#' @examples
#' removeConsecutiveDuplicates(c(1, 1, 2, 3, 3, 3, 4, 4, 5, 5))
#'
#'
removeConsecutiveDuplicates <- function(spreadValues){
  spreadValues <- as.matrix(spreadValues)
  uv <- unique(spreadValues)
  tra <- c()
  for(i in uv){
    ind <- which(spreadValues == i)
    tr <- which(diff(ind) == 1)
    if(length(tr)>0){
      tra <- c(tra,ind[tr])
    }
  }
  #
  if(length(tra) > 0){
    out <- spreadValues[-tra]
  }else{
    out <- spreadValues
  }
  return(out)
}



#' spreadPaleoData
#' @description A function to interpolate paleoData_values to a new vector using a nearest neighbor interpolation.
#'
#' @param age    age vector for original timeseries
#' @param value  paleoData vector to spread
#' @param newAge new age vector to spread over. If NA, this will be estimated using seq(min(age), max(age), spreadBy)
#' @param spreadBy if newAge is not supplied, what is the desired resolution of the spread values? If NA, the minimum age gap divided by 5 will be used
#' @param spreadMax a limit to how many years a value can be interpolated across (given in years). If NA, no limit is imposed.
#' @param spreadMaxPct a limit to how many years a value can be interpolated across (given in the quantile ranking of age gaps (0.5=the median age gap). A floor of 1 year is imposed). If NA, no limit is imposed. Superseded by maxGap
#' @param minAge minimum age value for output age vector. Default is -75 cal yr BP (2025 CE)
#' @param maxAge maximum age value for output age vector. If NA, this will be the max(age) plus half the gap between the two oldest ages
#'
#' @importFrom magrittr %>%
#'
#' @return list of length=2 (spreadAge and spreadVal numeric vectors)
#' @export
spreadPaleoData <- function(age,
                            value,
                            newAge = NA,
                            spreadBy = NA,
                            spreadMax = NA,
                            spreadMaxPct = NA,
                            minAge = -75,
                            maxAge = NA){

  #Check to make sure the input ages look ok
  if(length(age)==0){return(list(spreadAge = age,spreadVal = value))} #what's happening
  #remove nonfinite ages
  good <- which(is.finite(age))
  if(length(good) < 2){return(list(spreadAge = matrix(NA,ncol = length(age)),spreadVal = matrix(NA,ncol = length(age))))}

  #Remove any NA ages
  age   <- as.vector(age[good])
  value <- as.vector(value[good])
  #Sort ages consecutively
  ageidx <- order(age)
  age   <- age[ageidx]
  value <- value[ageidx]
  #if any duplicate ages, remove (calc average of values for same age)
  if (length(age) != (length(unique(age)))){
    agesnew <- sort(unique(age))
    vals <- c()
    for (a in agesnew){vals <- c(vals,mean(value[age==a],na.rm=T))}
    value <- vals
    age <- agesnew
  }

  #find max age to cutoff. This is an issue if the newAge vector is too long for the record
  maxAge <- max(age)+diff(age)[length(diff(age))]/2

  #get vector to spread data to
  if(any(is.na(newAge))){
    # if NA, use the minimum age gap divided by 5 to get a spreadBy value
    if (is.na(spreadBy)){spreadBy <- ceiling(min(abs(diff(age))/5))}
    # if NA, create a newAge which considered the age range of the proxy record
    # and adds additional buffer based on the resolution at either end
    newAge <- seq(ceiling(min(age)-diff(age)[1]/2),
                  floor(maxAge),
                  by = spreadBy)
  }else{
    spreadBy <- min(diff(sort(newAge)))
  }

  #Extend data as needed
  if(min(age) > min(newAge)){# we need to extend age
    age <- c(min(newAge),age)
    value <- c(value[1],value)
  }
  if(max(age) < max(newAge)){# we need to extend age
    age <- c(age,max(newAge))
    value <- c(value,value[length(value)])
  }

  #One-dimensional interpolation of points.
  newVals <- pracma::interp1(as.vector(age),as.vector(value),xi = newAge,method = 'nearest')

  #remove values that exceed maxGap
  if (!is.na(spreadMaxPct) & is.na(spreadMax)){
    spreadMax <- stats::quantile(diff(sort(age)),spreadMaxPct)
  }
  if(!is.na(spreadMax)){
    #distance to nearest
    d2n <- purrr::map_dbl(newAge,function(x) min(abs(x-age)))
    #assign NA where gap is too large
    newVals[d2n > spreadMax] <- NA
  }

  #remove values that are too young or old
  good <- dplyr::between(newAge,minAge,maxAge)
  if(sum(!good) > 0){
    newAge[!good] <- NA
    newVals[!good] <- NA
  }

  #return
  return(list(spreadAge = newAge,spreadVal = newVals))

}


#' simple binning of a lipd_ts object paleoData_values using a nearest neighbor approach
#'
#' @param ts a lipd_ts object
#' @param binvec vector of time boundaries over which to bin
#' @param ageVar specify the name the time variable (typically 'age' or 'year')
#' @param spread should values be interpolated between bins? (TRUE/FALSE)
#' @param gaussianizeInput Force values to gaussian distribution before analysis (TRUE/FALSE)
#' @param alignInterpDirection multiply values by -1 if the interpretation direction is negative (TRUE/FALSE). See [interpDirectionSign()].
#' @param scope the interpretation scope to align with (typically "climate" (default) or "isotope"). See [interpDirectionSign()].
#' @inheritParams spreadPaleoData
#'
#' @return numeric vector of values of equal length to binvec
#' @export
#' @importFrom pracma interp1
simpleBinTs <- function(ts,
                        binvec,
                        ageVar = "age",
                        spread = TRUE,
                        spreadBy = abs(mean(diff(binvec)))/10,
                        spreadMax = as.numeric(abs(stats::quantile(abs(diff(sort(ts[[ageVar]]))),probs = .75,na.rm = TRUE))),
                        gaussianizeInput= FALSE,
                        alignInterpDirection = TRUE,
                        scope = "climate"){

  #Check not ensemble

  if(NCOL(ts[[ageVar]])>1){stop('ages must be given to simpleBinTs() as a 1d vector not a matrix. Use sampleEnsembleThenBinTs() or calculate the ensemble median')}
  if(NCOL(ts[['paleoData_values']])>1){stop('paleoDataValues must be given to simpleBinTs() as a 1d vector not a matrix. Use sampleEnsembleThenBinTs() or calculate the ensemble median')}

  #Spread data if TRUE (default = TRUE)
  if(spread){#estimate for contiguous sampling with a nearest neighbor interpolation
    sp <- spreadPaleoData(age = ts[[ageVar]],
                          value = ts$paleoData_values,
                          spreadBy = spreadBy,
                          spreadMax = spreadMax)

    age <- sp$spreadAge
    vals <- sp$spreadVal
  }else{#use without any spreading
    age <- ts[[ageVar]]
    vals <-ts$paleoData_values
  }

  #gaussianize if TRUE (default = FALSE)
  if(gaussianizeInput){
    vals <- ens::gaussianize(vals)
  }

  #align if TRUE (default = TRUE)
  if(alignInterpDirection){
    di <- interpDirectionSign(ts, scope = scope)
    if(!is.na(di) && di < 0){
      vals <- vals * -1
    }
  }

  #bin
  binnedVals <- ens::bin(age,vals,bin.vec = binvec)[,2]

  return(binnedVals)
}



#' sampleEnsembleThenBinTs
#'
#' @inheritParams simpleBinTs
#' @param uncVar specify the name the uncertainty variable
#' @param defaultUnc a uncertainty to use if uncVar is NULL (default = 1.5, in the units given by `uncUnits`)
#' @param uncUnits units of `defaultUnc`: "native" (default) for the record's own units (paleoData_units), or "sd" for multiples of the record's standard deviation, so that `defaultUnc = 0.5` adds noise with half the record's SD whatever its units. Uncertainties read from the timeseries (`uncVar`) are measured values and are always in native units.
#' @param ar Autocorrelation coefficient to use for modelling uncertainty on paleoData, what fraction of the uncertainties are autocorrelated? (default = sqrt(0.5); or 50 percent autocorrelated uncertainty)
#' @param bamModel a list that describes the model to use in BAM (default = list(ns = 1, name = "bernoulli", param = 0.05))
#' @param exactAgeArchives archive types (matched case-insensitively against `ts$archiveType`) whose ages are treated as exact, so no age uncertainty is simulated for them when they have a single age vector. Defaults to tree rings ("Wood") and historical documents. Use `character(0)` to simulate age uncertainty for every record.
#'
#' @return  numeric vector of values of equal length to binvec
#' @export
sampleEnsembleThenBinTs <- function(ts,
                                    binvec,
                                    ageVar = "age",
                                    uncVar = "paleoData_uncertainty1sd",
                                    defaultUnc = 1.5,
                                    uncUnits = c("native", "sd"),
                                    ar = sqrt(0.5),
                                    bamModel = list(ns = 1, name = "bernoulli", param = 0.05),
                                    exactAgeArchives = c("Wood", "Documents"),
                                    spread = TRUE,
                                    spreadBy = abs(mean(diff(binvec)))/10,
                                    spreadMax =  as.numeric(stats::quantile(abs(diff(thisAge)),probs = .75,na.rm = TRUE)),
                                    gaussianizeInput = FALSE,
                                    alignInterpDirection = TRUE,
                                    scope = "climate"){
  ts_sampled <- ts
  #sample from ageEnsemble
  if(is.null(ts[[ageVar]])){
    stop(print(paste0(ts$dataSetName,": has a null for its age variable")))
  }
  if(NCOL(ts[[ageVar]]) > 1){#draw from ensemble
    thisAge <- ts[[ageVar]][ , sample.int(NCOL(ts[[ageVar]]),size = 1)]
  }else if(isExactlyDated(ts, exactAgeArchives)){#precisely dated (e.g. tree rings): use the ages as they are
    thisAge <- ts[[ageVar]]
  }else{#simulate
    thisAge <- ens::simulateBam(matrix(1,nrow = length(ts[[ageVar]])),as.matrix(ts[[ageVar]]),model = bamModel,ageEnsOut = TRUE)$ageEns
  }
  ts_sampled[[ageVar]] <- thisAge

  #Now sample from paleoData
  if(NCOL(ts$paleoData_values) > 1){#draw from ensemble
    thisPdv <- ts$paleoData_values[ , sample.int(NCOL(ts$paleoData_values),size = 1)]
  }else{ #simulate uncertainty from number
    if(!is.null(ts[[uncVar]])){#measured uncertainty, in the record's units
      tu <- as.numeric(ts[[uncVar]])
    }else if(match.arg(uncUnits) == "sd"){#default uncertainty, in standard deviations of the record
      tu <- defaultUnc * stats::sd(ts$paleoData_values, na.rm = TRUE)
    }else{
      tu <- defaultUnc
    }
    thisPdv <- ts$paleoData_values+simulateAutoCorrelatedUncertainty(sd = tu,n = length(ts$paleoData_values),ar = ar)
  }
  ts_sampled$paleoData_values <- thisPdv

  #check
  if(length(thisPdv) != length(thisAge)){stop("Paleodata and ages must have the same number of observations.")}

  #bin values using simpleBinTs()
  binnedVals <- simpleBinTs(ts_sampled, binvec, ageVar = ageVar, spread = spread,
                            spreadBy = spreadBy,
                            spreadMax = spreadMax,
                            gaussianizeInput = gaussianizeInput,
                            alignInterpDirection = alignInterpDirection,
                            scope = scope)

  return(binnedVals)

}


#' Sign of a timeseries' interpretation direction
#'
#' Finds the interpretation direction of a LiPD timeseries and returns 1 if
#' larger values mean more of the interpreted variable, -1 if they mean less,
#' and NA if there is no usable direction.
#'
#' Directions are read from `interpretationN_direction` fields (current LiPD),
#' or from legacy `climateInterpretationN_interpDirection` style fields.
#' Values beginning with "pos" or "incr" count as positive, "neg" or "decr"
#' as negative (case-insensitive), and numbers by their sign. Anything else
#' (e.g. "decrease/increase") gives NA.
#'
#' If `scope` is given, the first interpretation with that scope
#' (`interpretationN_scope`) is used. If none has that scope, the first
#' interpretation without a scope is used. If `scope` is NA, the first
#' interpretation with a direction is used.
#'
#' @param ts a single timeseries (one element of a TS list, or one row of a timeseries tibble)
#' @param scope interpretation scope to use, e.g. "climate" (default) or "isotope"
#'
#' @return 1, -1 or NA
#' @export
interpDirectionSign <- function(ts, scope = "climate"){
  if(is.data.frame(ts)){
    ts <- purrr::transpose(as.list(ts))[[1]]
  }
  nm <- names(ts)
  dirNames <- nm[grepl("^interpretation[0-9]+_direction$", nm) |
                   grepl("[Ii]nterpretation[0-9]+_interpDirection$", nm)]
  if(length(dirNames) == 0){
    return(NA_real_)
  }
  getVal <- function(n){
    v <- unlist(ts[[n]])
    if(length(v) == 0){
      return(NA)
    }
    v[1]
  }
  scopeOf <- function(dn){
    #legacy names carry the scope as a prefix, e.g. climateInterpretation1_interpDirection
    if(grepl("^[a-z]+Interpretation[0-9]+_interpDirection$", dn)){
      return(sub("Interpretation.*$", "", dn))
    }
    sn <- sub("_(interpD|d)irection$", "_scope", dn)
    s <- if(sn %in% nm) getVal(sn) else NA
    if(is.na(s) || !nzchar(trimws(as.character(s)))){
      return(NA_character_)
    }
    tolower(trimws(as.character(s)))
  }
  parseDir <- function(d){
    if(is.null(d) || length(d) == 0 || is.na(d)){
      return(NA_real_)
    }
    if(is.numeric(d)){
      return(if(d < 0) -1 else if(d > 0) 1 else NA_real_)
    }
    d <- tolower(trimws(as.character(d)))
    #ambiguous, e.g. "decrease/increase" or "t_air (positive), p_amount (negative)"
    if(grepl("[/,;]", d)){
      return(NA_real_)
    }
    if(grepl("^(pos|incr)", d)){
      return(1)
    }
    if(grepl("^(neg|decr)", d)){
      return(-1)
    }
    num <- suppressWarnings(as.numeric(d))
    if(!is.na(num) && num != 0){
      return(sign(num))
    }
    NA_real_
  }
  #order by interpretation number
  idx <- suppressWarnings(as.numeric(sub("^.*[Ii]nterpretation([0-9]+)_.*$", "\\1", dirNames)))
  dirNames <- dirNames[order(idx)]
  dirs <- vapply(dirNames, function(n) parseDir(getVal(n)), numeric(1))
  if(is.na(scope)){
    ok <- which(!is.na(dirs))
    return(if(length(ok)) unname(dirs[ok[1]]) else NA_real_)
  }
  scopes <- vapply(dirNames, scopeOf, character(1))
  hit <- which(scopes == tolower(scope))
  if(length(hit) == 0){
    hit <- which(is.na(scopes))
  }
  if(length(hit) == 0){
    return(NA_real_)
  }
  unname(dirs[hit[1]])
}

#' Is a timeseries from an archive with exact ages?
#' @param ts a single timeseries
#' @param exactAgeArchives archive types treated as exactly dated
#' @return TRUE/FALSE
#' @keywords internal
isExactlyDated <- function(ts, exactAgeArchives){
  if(length(exactAgeArchives) == 0){
    return(FALSE)
  }
  at <- unlist(ts$archiveType)
  if(length(at) == 0 || is.na(at[1])){
    return(FALSE)
  }
  tolower(at[1]) %in% tolower(exactAgeArchives)
}
