library(lipdR) #remotes::install_github("nickmckay/lipdR)
library(tidyverse)
library(purrr)
library(geoChronR) #remotes::install_github("nickmckay/geoChronR)

#load the database
D <- readLipd("https://lipdverse.org/Temp12k/1_0_2/Temp12k1_0_2.zip")

# create a tibble to work with
tts <- ts2tibble(sTS)


#filter the database to get just the data you want. At the very least you probably want to get just the Temp12k series
temp12k <- filter(tts,paleoData_inCompilation == "Temp12k")

#Decide how you will get the data onto a common timescale
bin.vec <- seq(-50,12000,by = 200)
bin.mid <- data.frame(binMiddle = rowMeans(cbind(bin.vec[-1],bin.vec[-length(bin.vec)])))
binFun <- function(age,vals,...){geoChronR::bin(age,vals,...)$y}

#Apply that function to all the selected datasets
dataMatrix <- map2_dfc(temp12k$age,temp12k$paleoData_values,binFun,bin.vec) %>%
  setNames(paste0(temp12k$dataSetName,"-",temp12k$paleoData_variableName))

#Then add time as the first column
final <- bind_cols(bin.mid,dataMatrix,.name_repair = "minimal")





