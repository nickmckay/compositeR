standardizeOverInterval <- function(ages,pdm,interval,minN = 8,normalizeVariance = TRUE){

  #check the data
  if(length(ages)!=nrow(pdm)){
    stop("the ages must match the rows in the paleoData matrix")
  }

  #get the interval
  iStart <- min(interval)
  iEnd <- max(interval)

  #check which cols have enough data in
  wa <- which(ages>=iStart & ages<iEnd)
  if(sum(is.finite(ages[wa])))

  #subset the matrix by the interval


  spdm <- pdm[wa,]
  sm <- apply(spdm,2,mean,na.rm = TRUE)
  ss <- apply(spdm,2,sd,na.rm = TRUE)


  #check which cols have enough data in the window

  ni <- apply(spdm,2,function(x) sum(is.finite(removeConsecutiveDuplicates(x))))

  badc <- which(ni<minN)

 pdm[,badc] <- NA
 spdm[,badc] <- NA

  #scale the matrix
  if(normalizeVariance){
    scaledPdm <- scale(pdm,center = sm,scale = ss)
  }else{
    scaledPdm <- scale(pdm,center = sm,scale = FALSE)
  }


  #return it
  return(scaledPdm)
}

standardizeOverRandomInterval <- function(ages,pdm,duration,searchRange,normalizeVariance = TRUE,minN = 8){

  #check the data
  if(length(ages)!=nrow(pdm)){
    stop("the ages must match the rows in the paleoData matrix")
  }

  scaledPdm <- matrix(NA,nrow = nrow(pdm),ncol = ncol(pdm))
  for(i in 1:ncol(pdm)){

    #get the possible range

    goodAges <-   ages[!is.na(pdm[,i])]

    pStart <- min(goodAges)
    pEnd <- max(goodAges)

    if(pEnd-pStart < duration){
      warning(paste("column",i,"doesnt have the required duration"))
      scaledPdm[,i] <- NA
      next
    }

    startOptions <- goodAges[goodAges>= max(min(c(searchRange,pStart))) & goodAges<=min(c(max(searchRange),pEnd)-duration)]

    if(length(startOptions) < 1){
      warning(paste("column",i,"doesnt have an overlap between the required duration and the searchRange"))
      scaledPdm[,i] <- NA
      next
    }

    nVal <- 0
    tt <- 0
    posStarts <- sample(startOptions)
    for(ps in posStarts){
      iStart <- ps
      iEnd <- iStart+duration

      #check nvalues
      pass <- which(ages>=iStart & ages<iEnd)
      nVal <- sum(is.finite(removeConsecutiveDuplicates(pdm[pass,i])))

      if(nVal >= minN){
        break #move on
      }
    }

    if(nVal < minN){
      warning(paste("cant find an interval in column",i,"that has enough (at least minN (",minN,") observations"))
      scaledPdm[,i] <- NA
      next
    }

    #subset the matrix by the interval
    spdm <- as.matrix(pdm[pass,i])

    sm <- mean(spdm,na.rm = TRUE)
    ss <- sd(removeConsecutiveDuplicates(spdm),na.rm = TRUE)
    #scale the matrix
    if(normalizeVariance){
      scaledPdm[,i] <- scale(pdm[,i],center = sm,scale = ss)
    }else{
      scaledPdm[,i] <- scale(pdm[,i],center = sm,scale = FALSE)
    }

  }

  #return it
  return(scaledPdm)
}







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

standardizeMeanIteratively <- function(ages,pdm,duration,searchRange,normalizeVariance = TRUE,thresh = 0.01,minN = 8){
  #nicks crackpot idea

  #start by scaling everything

  start <- standardizeOverRandomInterval(ages = ages,pdm = pdm,duration = duration,searchRange = searchRange,normalizeVariance = normalizeVariance,minN = minN)

  #remove records that failed (this should potentially cause an error in the future)
  start[!is.finite(start)] <- NA
  colsds <- apply(pdm,2,sd,na.rm =TRUE)
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






