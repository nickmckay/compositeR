## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

getCompilationVersions <- function(ts,compilation,version){
  out <- FALSE
  for (comp in names(ts)[grepl('compilationName',names(ts),ignore.case = T)]){
    name <- ts[[comp]]
    if (tolower(name) == tolower(compilation)){
      versions <- ts[[stringr::str_replace(comp, 'Name', 'Version')]]
      if (version %in% versions){
        out <- TRUE
      }
    }
  }
  return(out)
}


version <- '1_0_2'
compilation <- 'Temp12k'
#Load LiPD Data
D <- lipdR::readLipd(file.path("https://lipdverse.org",compilation,version,paste0(compilation,version,'.zip')))
#Extract timeseries objects
mts <- lipdR::splitInterpretationByScope(lipdR::extractTs(D))
#Filter for temp12k
Temp12k <- mts[unlist(purrr::map(mts,getCompilationVersions,compilation,version))]
#Save filtered compilation
usethis::use_data(Temp12k)


version <- '1_0_0'
compilation <- 'temp12kEnsemble'
#Load LiPD files with age ensembles
#from https://github.com/nickmckay/Temperature12k/tree/ea854570be681680dc724fd6b2f31c0f4667483c/ScientificDataAnalysis/lipdFilesWithEnsembles
lipdFilesWithEnsembles <- lipdR::readLipd(path="/data-raw/lipdFilesWithEnsembles/")
#usethis::use_data(lipdFilesWithEnsembles)

#Extract timeseries objects
mts <- lipdR::splitInterpretationByScope(lipdR::extractTs(lipdFilesWithEnsembles))
#Filter for temp12k
Temp12k_ageens <- mts[which(lipdR::pullTsVariable(mts,"paleoData_inCompilation")=="temp12kEnsemble")]#unlist(purrr::map(mts,getCompilationVersions,compilation,version))]
Temp12k_ageens <- Temp12k_ageens[which(lipdR::pullTsVariable(Temp12k_ageens,'climateInterpretation1_seasonalityGeneral') %in% c("annual","winterOnly","summerOnly"))]

#Save filtered compilation
usethis::use_data(Temp12k_ageens)

