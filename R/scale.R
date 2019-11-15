scaleComposite <- function(composite,binvec,scaleYears,scaleData,scaleWindow = NA,rescale = TRUE){

  if(NCOL(scaleData) > 1){#ensemble! :)
    d <- bin(scaleYears,values = scaleData[,sample.int(ncol(scaleData),size = 1)],binvec = binvec)

  }else{ #nonsemble :(
    d <- bin(scaleYears,values = scaleData,binvec = binvec)

  }

  #restrict the window (Use the whole scale data if not)
  if(is.na(scaleWindow)){
    scaleWindow <- range(scaleYears)
  }
    good <- which(d$x >= min(scaleWindow) & d$x <= max(scaleWindow))
    dv <- d$y[good]

  m <- mean(dv,na.rm = TRUE)
  if(rescale){#forces the mean of the final scaled window to be 0
    m <- 0
  }


  s <- sd(dv,na.rm = TRUE)

  #now scale the composite to mean = 0 and sd = 1 over scale Window
  compYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))
  swp <- which(compYears >= min(scaleWindow) & compYears <= max(scaleWindow))
  scp <- scale(composite,center = mean(composite[swp],na.rm = TRUE),scale  = sd(composite[swp],na.rm = TRUE))

  #rescale
  scaled <- as.matrix(scp)*s+m

return(scaled)

}

