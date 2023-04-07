#' remove consecutive duplicates
#'
#' @param spreadValues
#'
#' @return
#' @export
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


  if(length(tra) > 0){
    out <- spreadValues[-tra]
  }else{
    out <- spreadValues
  }
  return(out)


}




#' spreadPaleoData
#'
#' @param age
#' @param value
#' @param newAge new age vector to spread over
#' @param spreadBy if newAge is not supplied, what is the desired resolution of spread data.
#' @param maxGap
#' @param maxPct
#' @param minAge
#'
#' @return
#' @export
spreadPaleoData <- function(age,
                            value,
                            newAge = NA,
                            spreadBy,
                            maxGap = NA,
                            maxPct = 0.75,
                            minAge = -69){

  if(length(age)==0){
    #what's happening
    return(list(spreadAge = age,spreadVal = value))
  }

  #remove nonfinite ages
  good <- which(is.finite(age))
  if(length(good) < 2){
    return(list(spreadAge = matrix(NA,ncol = length(age)),spreadVal = matrix(NA,ncol = length(age))))
  }
  age <- age[good]
  value <- value[good]


  hasNas <- FALSE
  # #it doesn't handle NAs appropriately, so lets fix that
  # if(any(is.na(value))){
  #   hasNas <- TRUE
  #   value[is.na(value)] <- -999999
  # }

  #spread and interpolate
  if(all(is.na(newAge))){
    newAge <- seq(ceiling(min(age)),floor(max(age)),by = spreadBy)
  }else{
    spreadBy <- median(newAge,na.rm = TRUE)
  }


  age <- as.vector(age)
  as <- sort(age,index.return = TRUE)
  age <- age[as$ix]
  value <- as.vector(value)
  value <- value[as$ix]
  if(min(age) > min(newAge)){# we need to extend age
    age <- c(min(newAge),age)
    value <- c(value[1],value)
  }

  if(max(age) < max(newAge)){# we need to extend age
    age <- c(age,max(newAge))
    value <- c(value,value[length(value)])
  }

  newVals <- pracma::interp1(as.vector(age),as.vector(value),xi = newAge,method = "nearest")

  # if(hasNas){
  #   newVals[which(newVals == -999999)] <- NA
  # }


  #add on to the beginning
  nv1 <- which(newVals == value[1])
  if(any(diff(nv1) != 1)){
    nv1 <- nv1[1:(min(which(diff(nv1) != 1)))]
  }


  f <- min(newAge)-min(c(0,length(nv1)-1))*spreadBy
  t <- min(newAge)-spreadBy

  if(f < t){
  begAge <- seq(from = f,to = t,by = spreadBy)
  begVal <- rep(newVals[1],times = length(begAge))
  }else{
    begAge <- c()
    begVal <- c()
  }

  #add on to the end
  end <- length(newAge)
  nv2 <- which(newVals == value[length(value)])
  if(any(diff(nv2) != 1)){
    nv2 <- nv2[max(which(diff(nv2) != 1)):length(nv2)]
  }


  f <- max(newAge)+spreadBy
  t <- max(newAge)+max(c(0,length(nv2)-1))*spreadBy
  if(f < t){
    endAge <- seq(from = f,to = t,by = spreadBy)
    endVal <- rep(newVals[length(newVals)],times = length(endAge))
  }else{
    endAge <- c()
    endVal <- c()
  }



  #append them
  newAgeOut <- c(begAge,newAge,endAge)
  newValsOut <- c(begVal,newVals,endVal)


  #distance to nearest
  d2n <- map_dbl(newAgeOut,function(x) min(abs(x-age)))

  #get local d2n maxima
  locmaxi <- which(diff(sign(diff(d2n)))==-2)+1
  if (length(locmaxi)>0){
    locmax <- pracma::interp1(x = as.vector(c(newAgeOut[1],newAgeOut[locmaxi],newAgeOut[length(newAgeOut)])),
                              as.vector(c(d2n[locmaxi[1]],d2n[locmaxi],d2n[locmaxi[length(locmaxi)]]))
                              ,xi = newAgeOut,
                              method = "nearest")
    lmpct <- d2n/locmax
  } else{
    lmpct <-0
  }


  #remove values that exceed maxGap
  if(!is.na(maxGap)){
    newValsOut[d2n > maxGap] <- NA
  }

  if(!is.na(maxPct)){
    newValsOut[lmpct > maxPct] <- NA
  }

  #remove values that are too young
  ty <- which(newAgeOut < minAge)
  if(length(ty) > 0){
    newAgeOut <- newAgeOut[-ty]
    newValsOut <- newValsOut[-ty]
  }


  return(list(spreadAge = newAgeOut,spreadVal = newValsOut))

}


#' simple binning of a TS object instance
#'
#' @param ts a ts object
#' @param binvec vector of boundaries
#' @param ageVar age variable
#'
#' @return
#' @export
#' @importFrom pracma interp1
#'
#' @examples
simpleBinTs <- function(ts,
                        binvec,
                        ageVar = "age",
                        spread = TRUE,
                        spreadBy = abs(mean(diff(binvec)))/10,
                        gaussianizeInput= FALSE,
                        alignInterpDirection = TRUE,
                        scope = "climate"){
  if(spread){#estimate for contiguous sampling with a nearest neighbor interpolation

  sp <- spreadPaleoData(age = ts[[ageVar]],
                        value = ts$paleoData_values,
                        spreadBy = spreadBy,
                        maxGap = as.numeric(quantile(abs(diff(ts[[ageVar]])),probs = .75,na.rm = TRUE)))
  age <- sp$spreadAge
  vals <- sp$spreadVal


  }else{#use without any spreading
    age <- ts[[ageVar]]
    vals <-ts$paleoData_values
  }

  #gaussianize?
  if(gaussianizeInput){
    vals <- geoChronR::gaussianize(vals)
  }


  if(alignInterpDirection){
    #check for direction
    din <- names(ts)[stringr::str_detect("_interpDirection",string = names(ts))]
    di <- unlist(magrittr::extract(ts,din))

    sin <- names(ts)[stringr::str_detect("_scope",string = names(ts))]
    si <- unlist(magrittr::extract(ts,sin))

    if(!is.na(scope)){
      di <- di[grepl(pattern = scope,x = si)]
    }

    if(length(di)>0){
      if(all(grepl(di,pattern = "negative",ignore.case = TRUE))){
        vals <- vals * -1
      }
    }
  }

  bd <- geoChronR::bin(age,vals,bin.vec = binvec)[,2]
  return(bd)
}



#' sampleEnsembleThenBinTs
#'
#' @param ts
#' @param binvec
#' @param ageVar
#' @param uncVar
#' @param defaultUnc
#' @param ar
#' @param bamModel
#' @param spread
#' @param spreadBy
#' @param gaussianizeInput
#' @param alignInterpDirection
#' @param scope
#'
#' @return
#' @export
sampleEnsembleThenBinTs <- function(ts,
                                    binvec,
                                    ageVar = "age",
                                    uncVar = "paleoData_uncertainty1sd",
                                    defaultUnc = 1.5,
                                    ar = sqrt(0.5),
                                    bamModel = list(ns = 1, name = "bernoulli", param = 0.05),
                                    spread = TRUE,
                                    spreadBy = abs(mean(diff(binvec)))/10,
                                    gaussianizeInput = FALSE,
                                    alignInterpDirection = TRUE,
                                    scope = "climate"){
  #sample from ageEnsemble
  if(is.null(ts[[ageVar]])){
    stop(print(paste0(ts$dataSetName,": has a null for its age variable")))
  }
  if(NCOL(ts[[ageVar]]) > 1){#draw from ensemble
    thisAge <- ts[[ageVar]][ , sample.int(NCOL(ts[[ageVar]]),size = 1)]
  }else{#simulate
    thisAge <- geoChronR::simulateBam(matrix(1,nrow = length(ts[[ageVar]])),as.matrix(ts[[ageVar]]),model = bamModel,ageEnsOut = TRUE)$ageEns
  }


  #Now sample from paleoData
  if(NCOL(ts$paleoData_values) > 1){#draw from ensemble
    thisPdv <- ts$paleoData_values[ , sample.int(NCOL(ts$paleoData_values),size = 1)]
  }else{ #simulate uncertainty from number
    if(!is.null(ts[[uncVar]])){
      tu <- ts[[uncVar]]
    }else{
      tu <- defaultUnc
    }
    thisPdv <- ts$paleoData_values+simulateAutoCorrelatedUncertainty(sd = tu,n = length(ts$paleoData_values),ar = ar)
  }

  #check
  if(length(thisPdv) != length(thisAge)){
    stop("Paleodata and ages must have the same number of observations.")
  }

  #gaussianize?
  if(gaussianizeInput){
    thisPdv <- geoChronR::gaussianize(thisPdv)
  }


  #spread?
  if(spread){#estimate for contiguous sampling with a nearest neighbor interpolation

    sp <- spreadPaleoData(age = thisAge,
                          value = thisPdv,
                          spreadBy = spreadBy,
                          maxGap = as.numeric(quantile(abs(diff(thisAge)),probs = .75,na.rm = TRUE)))
    age <- sp$spreadAge
    vals <- sp$spreadVal

  }else{#use without any spreading
    age <- thisAge
    vals <-thisPdv
  }


  if(alignInterpDirection){
    #check for direction
    din <- names(ts)[stringr::str_detect("_interpDirection",string = names(ts))]
    di <- unlist(magrittr::extract(ts,din))

    sin <- names(ts)[stringr::str_detect("_scope",string = names(ts))]
    si <- unlist(magrittr::extract(ts,sin))

    if(!is.na(scope)){
    di <- di[grepl(pattern = scope,x = si)]
    }

    if(length(di)>0){
      if(all(grepl(di,pattern = "negative",ignore.case = TRUE))){
        vals <- vals * -1
      }
    }
  }

  bd <- geoChronR::bin(time = age,values = vals,bin.vec = binvec)[,2]
  return(bd)

}
