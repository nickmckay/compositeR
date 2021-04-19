library(geoChronR)
library(lipdR)
library(purrr)
library(magrittr)
#load database
#load("~/GitHub/lipdverse/html/iso2k/current_version/iso2k0_14")

 D <- readLipd("https://lipdverse.org/iso2k/current_version/iso2k1_0_0.zip")

# #extract timeseries
TS <- extractTs(D)
sTS <- splitInterpretationByScope(TS)

tibTS <- ts2tibble(TS)

#get all the names in the TS
allnames <- sapply(sTS,names) %>%
  unlist() %>%
  unique() %>%
  sort()

names(tibTS) %>% sort()


#pull variables needed for filtering
primary <- pullTsVariable(sTS,variable = "paleoData_iso2kPrimaryTimeseries")
img <- pullTsVariable(sTS,variable = "paleoData_inferredMaterialGroup")
ivg <- pullTsVariable(sTS,variable = "isotopeInterpretation1_variableGroup")


#do filtering
temp <- which(ivg == "Temperature")

fTS <- sTS[intersect(temp,which(primary == "TRUE"))]


ls <- map_dbl(fTS,function(x) sum(!is.na(x$paleoData_values) & !is.na(x$year)))
ls2 <- map_dbl(fTS,function(x) length(x$paleoData_values))

fTS <- fTS[which(ls > 10 & ls2 >10)]



#bin the TS
binvec <-  seq(1, to = 2000)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))


#setup ensemble
nens <- 20

library(foreach)
library(doParallel)
registerDoParallel(2)

ensOut <- foreach(i = 1:nens) %dopar% {
    tc <- compositeEnsembles(fTS,
                             binvec,
                             stanFun = standardizeMeanIteratively,
                             binFun = simpleBinTs,
                             ageVar = "year",
                             alignInterpDirection = FALSE,
                             spread = TRUE,
                             duration = 50,
                             searchRange = c(0,2000),
                             normalizeVariance = FALSE)

  return(list(composite = tc$composite,count = tc$count))
}

# #plotting!

thisComposite <-  as.matrix(purrr::map_dfc(ensOut,extract,"composite"))
library(ggplot2)
compPlot <- plotTimeseriesEnsRibbons(X = binYears,Y = thisComposite)+
  scale_x_continuous(name = "Year (AD)",oob = scales::squish)+
  scale_y_continuous(name = "d18O (permil)",oob = scales::squish)+
  theme_bw()+
  ggtitle("Temperature sensitive")

compPlot
