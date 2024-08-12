#' standardizeOverInterval
#'
#' @param ages vector of ages. Each age corresponds to a row in pdm
#' @param pdm matrix of binned paleoData values
#' @param interval the age range to use for standardizing the data
#' @param minN minimum number length of de within the interval to calculate
#' @param normalizeVariance Should the pdm be scaled to a uniform mean (FALSE) or variance and mean (TRUE) (default = TRUE)
#'
#'
#' @return a matrix of standardized values
#' @export
standardizeOverInterval <- function(ages,pdm,interval,minN = 8,normalizeVariance = TRUE){

  #check the data
  if(length(ages)!=nrow(pdm)){stop("the ages must match the rows in the paleoData matrix")}

  #check which cols have enough data in
  wa <- which(ages>=min(interval) & ages<max(interval))
  if(sum(is.finite(ages[wa])) == 0){
    warning("binvec ages are outside of standarization interval")
    pdm[,] <- NA
    return(pdm)
  }

  #subset the matrix by the interval
  spdm <- as.matrix(pdm[wa,])

  #check which cols have enough data in the window (and assing as NA if )
  ni <- apply(spdm,2,function(x) sum(is.finite(removeConsecutiveDuplicates(x))))
  pdm[,which(ni<minN)] <- NA
  spdm[,which(ni<minN)] <- NA

  #Get mean and standard deviation for interval
  sm <- apply(spdm,2,mean,na.rm = TRUE)
  ss <- apply(spdm,2,stats::sd,na.rm = TRUE)

  #scale the matrix
  if(normalizeVariance){
    scaledPdm <- base::scale(pdm,center = sm,scale = ss)
  }else{
    scaledPdm <- base::scale(pdm,center = sm,scale = FALSE)
  }

  #return it
  return(scaledPdm)
}


#' standardizeOverRandomInterval
#' @param searchRange the age range within which to search for an interval
#' @param duration length of the interval to standardize the data within. Must be <= than searchRange
#' @inheritParams standardizeOverInterval
#'
#' @return a matrix of standardized values
#' @export
standardizeOverRandomInterval <- function(ages,pdm,duration,searchRange,minN = 8,normalizeVariance = TRUE){

  if (duration>diff(searchRange)){stop("The duration must be <= the age range indicated by searchRange")}

  #Find a random searchRange to use
  interval <- sample(seq(min(searchRange),max(searchRange)-duration,1),1)
  interval <- c(interval,interval+duration)

  scaledPdm <- standardizeOverInterval(ages, pdm, interval=interval, minN = minN, normalizeVariance = normalizeVariance)

  # Simplified this code (for better or worse) so that the same interval is chosen. If a record fails minN, so be it -Chris

  # #check the data
  # if(length(ages)!=nrow(pdm)){
  #   stop("the ages must match the rows in the paleoData matrix")
  # }
  # scaledPdm <- matrix(NA,nrow = nrow(pdm),ncol = ncol(pdm))
  # for(i in 1:ncol(pdm)){
  #
  #   #get the possible range
  #
  #   goodAges <-   ages[!is.na(pdm[,i])]
  #
  #   pStart <- min(goodAges)
  #   pEnd <- max(goodAges)
  #
  #   if(pEnd-pStart < duration){
  #     warning(paste("column",i,"doesnt have the required duration"))
  #     scaledPdm[,i] <- NA
  #     next
  #   }
  #
  #   startOptions <- goodAges[goodAges>= max(min(c(searchRange,pStart))) & goodAges<=min(c(max(searchRange),pEnd)-duration)]
  #
  #   if(length(startOptions) < 1){
  #     warning(paste("column",i,"doesnt have an overlap between the required duration and the searchRange"))
  #     scaledPdm[,i] <- NA
  #     next
  #   }
  #
  #   nVal <- 0
  #   tt <- 0
  #   posStarts <- sample(startOptions)
  #   for(ps in posStarts){
  #     iStart <- ps
  #     iEnd <- iStart+duration
  #
  #     #check nvalues
  #     pass <- which(ages>=iStart & ages<iEnd)
  #     nVal <- sum(is.finite(removeConsecutiveDuplicates(pdm[pass,i])))
  #
  #     if(nVal >= minN){
  #       break #move on
  #     }
  #   }
  #
  #   if(nVal < minN){
  #     warning(paste("cant find an interval in column",i,"that has enough (at least minN (",minN,") observations"))
  #     scaledPdm[,i] <- NA
  #     next
  #   }
  #
  #   #subset the matrix by the interval
  #   spdm <- as.matrix(pdm[pass,i])
  #
  #   sm <- mean(spdm,na.rm = TRUE)
  #   ss <- stats::sd(removeConsecutiveDuplicates(spdm),na.rm = TRUE)
  #   #scale the matrix
  #   if(normalizeVariance){
  #     scaledPdm[,i] <- scale(pdm[,i],center = sm,scale = ss)
  #   }else{
  #     scaledPdm[,i] <- scale(pdm[,i],center = sm,scale = FALSE)
  #   }
  #
  # }

  #return it
  return(scaledPdm)
}





# I'm too confused by the below at the moment to try and change it -CH. Defaulting to revised standardizeOverRandomInterval()

#' recordRMSE
#'
#' @param td output of standardizeOverRandomInterval
#' @param palMat matrix of paleoData values
#'
#' @return list with RMSE and Bias vectores
#' @export
recordRMSE <- function(td,palMat){
  #difference matrix
  tdm <- matrix(td,nrow = nrow(palMat),ncol = ncol(palMat),byrow = FALSE)
  dm <- tdm-palMat

  #remove identical columns
  good <- which(abs(colMeans(dm,na.rm = TRUE)) > .Machine$double.eps )

  #calculate row bias
  gm <- dm[,good]
  gm[!is.finite(gm)] <- NA

  #rowMeanBias <- rowMeans(gm,na.rm = TRUE)

  #calculate total bias
  totalBias <- mean(c(gm),na.rm = TRUE)

  #calculate total RMSE (of all values, not rows)
  totalRMSE <- sqrt(mean(c(gm^2),na.rm = TRUE))

  return(list(totalRMSE = totalRMSE, totalBias = totalBias))

}

optimizeMeanAndSd <- function(msd,td,palMat){
  m <- msd[1]
  s <- msd[2]

  ntd <- scale(td,center = m,scale = s)

  return(recordRMSE(ntd,palMat)$totalRMSE)
}

optimizeMean <- function(m,td,palMat){

  ntd <- scale(td,center = m,scale = FALSE)

  return(recordRMSE(ntd,palMat)$totalRMSE)
}


#' standardizeMeanIteratively
#'
#' @param ages vector of ages. Each age corresponds to a row in pdm
#' @param pdm matrix of binned paleoData values
#' @param duration length of the interval to standardize the data within
#' @param searchRange the age range to which the interval range will be sampled from
#' @param normalizeVariance Should the pdm be scaled to a uniform mean (FALSE) or variance and mean (TRUE) (default = TRUE)
#' @param thresh threshold for mean RMSE
#' @param minN minimum number length of de within the interval to calculate
#'
#' @return a matrix of standardized values
#' @export
standardizeMeanIteratively <- function(ages,
                                       pdm,
                                       duration,
                                       searchRange,
                                       normalizeVariance = TRUE,
                                       thresh = 0.01,
                                       minN = 8){
  #nicks crackpot idea

  #start by scaling everything

  start <- standardizeOverRandomInterval(ages = ages,
                                         pdm = pdm,
                                         duration = duration,
                                         searchRange = searchRange,
                                         normalizeVariance = normalizeVariance,
                                         minN = minN)

  #remove records that failed (this should potentially cause an error in the future)
  start[!is.finite(start)] <- NA
  colsds <- apply(pdm,2,stats::sd,na.rm =TRUE)
  filledBins <- apply(pdm,2,function(x) sum(is.finite(x)))

  bad <- which(colSums(is.finite(start))==0 | colsds < 0.1 | filledBins < minN)

  if(length(bad) >= (NCOL(start)-2)){
  stop("No good columns after standardization")
  }



  #set up while loop
  delta <- thresh+1
  meanAllRMSE <- 1000 #pick something large
  upmat <- start

  while(delta > thresh){
    allRMSE <- matrix(NA,nrow = ncol(start))

    rvec <- sample(1:ncol(start)) #randomize the order

    for(j in rvec){
      #optimize
      #r <- optim(par = mean(upmat[,j],na.rm = TRUE),fn = optimizeMean,td = upmat[,j],palMat = upmat)
      es <- recordRMSE(upmat[,j],upmat)
      meanBias <- es$totalBias
      #update values
      upmat[,j] <- upmat[,j]-meanBias
      allRMSE[j] <- es$totalRMSE
    }
    old <- meanAllRMSE
    meanAllRMSE <- mean(allRMSE,na.rm = TRUE)
    delta <- old-meanAllRMSE
  }

  return(upmat)
}






