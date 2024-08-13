devtools::load_all()
library(lipdR)
library(magrittr)
library(geoChronR)
library(tidyverse)
library(sf)
library(dggridR)

#Load data
project <- 'RAW'
dir <- '/Users/chrishancock/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/Research/Manuscript/'
RAWts     <- readRDS(file.path(dir,project,'Data','RAW_LiPDts.RDS'))
shiftsAll <- readRDS(file.path(dir,project,'Data',paste(project,'meanShifts.RDS',sep='_')))

#Filter data
shifts <- shiftsAll %>% filter(paleoData_units=='degC') %>% filter(interpVar %in% c('Temperature')) %>% filter(interpSeason %in% c('Annual','Summer','Winter'))
fTS <- RAWts[(pullTsVariable(RAWts,'paleoData_TSid') %in% shifts$paleoData_TSid)]
binvec <- sort(unique(c(shifts$time_start,shifts$time_end)))
binyrs <- sort(unique(c(shifts$time_mid)))

#find weights
dggs <- dgconstruct(res=4,pole_lat_deg=90,topology='HEXAGON')
polygons <- dggridR::dgearthgrid(dggs)
#locate each record's cell
shifts <- shifts %>% dplyr::mutate(cell=(dggridR::dgGEO_to_SEQNUM(dggs,geo_longitude,geo_latitude)$seqnum))
#assign regions
defineRegion <- function(latmin,latmax,lonmin,lonmax){
  nlon = length(seq(lonmin,lonmax,1))
  nlat = length(seq(latmin,latmax,1))
  Poly_Coord_df <- data.frame(lon= c(seq(lonmin,lonmax,1),rev(seq(lonmin,lonmax,1))),
                              lat= c(rep(latmin,nlon),rep(latmax,nlon))
  ) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    summarise((geometry = sf::st_combine(geometry))) %>%
    sf::st_cast("POLYGON")
  return(Poly_Coord_df)
}
regions <- rbind(defineRegion(58,90,-180,-100),
                 defineRegion(58,90,-100,-10),
                 defineRegion(58,71,-10,50),
                 defineRegion(71,90,-10,50),
                 defineRegion(58,90,50,110),
                 defineRegion(58,90,110,180))
sf::st_geometry(regions) <- "geometry"
regions$name <- c("North America","Greenland","Scandinavia","Svalbard","Central Russia","Eastern Russia")
shifts$region <- sf::st_join(st_as_sf(x =data.frame(lat=shifts$geo_latitude,lon=shifts$geo_longitude), coords = c("lon", "lat"), crs = st_crs(regions)),regions)$name


#add weights which consider number of record within each cell and the number of records at each site
shifts$weight <- NA
for (c in unique(shifts$cell)){
  cellShifts <- (shifts %>% filter (cell == c))
  cellweight <- 1/length(unique(cellShifts$dataSetName))
  for (s in unique(cellShifts$site)){
    proxyweight <- cellweight/length(unique((cellShifts %>% filter(site == s))$paleoData_TSid))
    shifts$weight[which((shifts$cell==c) & (shifts$site==s))] <- proxyweight
  }
}

# Add weight and climateInterpretation1_interpDirection (which is missing (or labeled differently for some records))
for (i in 1:length(fTS)){
  sel <- (shifts%>%filter(paleoData_TSid==fTS[[i]]$paleoData_TSid))
  fTS[[i]]$climateInterpretation1_interpDirection <- sel$interpDir[1]
  fTS[[i]]$geo_weight <- sel$weight[1]
  fTS[[i]]$geo_cell <- sel$cell[1]
  fTS[[i]]$geo_region <- sel$region[1]
}

ensOut_all <- list()
for (reg in c(regions$name,'all')){
  #do composite
  if (reg == 'all'){
    fTS_sel <- fTS
  }else{
    fTS_sel <- fTS[pullTsVariable(fTS,'geo_region')==reg]
  }
  ensOut <- compositeEnsembles2(
    fTS = fTS_sel,
    nens = 100,
    binvec = binvec,
    binFun = sampleEnsembleThenBinTs,
    ageVar = "ageEnsemble",
    uncVar = "paleoData_temperature12kUncertainty",
    weights = "geo_weight",
    stanFun = standardizeOverRandomInterval,
    searchRange = c(1000,11000),
    duration = 4000,
    normalizeVariance = TRUE,
    gaussianizeInput = TRUE,
    scale = FALSE,
    spread = TRUE,
    samplePct = 1,
  )
  print(reg)
  plot(ensOut)
  ensOut_all[[reg]] <- ensOut
}
#print(ensOut)


ensOut_allCal <- list()
for (reg in c('all')){
  #do composite
  print(reg)

  if (reg == 'all'){
    fTS_sel <- fTS[(pullTsVariable(fTS,'paleoData_units')=="degC")]
  }else{
    fTS_sel <- fTS[(pullTsVariable(fTS,'geo_region')==reg) & (pullTsVariable(fTS,'paleoData_units')=="degC")]
  }
  ensOut <- compositeEnsembles2(
    fTS = fTS_sel,
    nens = 100,
    binvec = binvec,
    binFun = sampleEnsembleThenBinTs,
    ageVar = "ageEnsemble",
    uncVar = "paleoData_temperature12kUncertainty",
    weights = "geo_weight",
    stanFun = standardizeOverRandomInterval,
    searchRange = c(0,12000),
    duration = 6000,
    normalizeVariance = FALSE,
    gaussianizeInput = FALSE,
    scale = FALSE,
    spread = TRUE,
    samplePct = 1,
  )
  plot(ensOut)
  ensOut_allCal[[reg]] <- ensOut
}
saveRDS(ensOut_allCal,file='/Users/chrishancock/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/Research/Manuscript/RAW/Data/proxyComposite_AnnualPlus')
#print(ensOut)

plotlist <- list()
for (i in (names(ensOut_allCal))){
  plotlist[[i]] <- plot(ensOut_allCal[[i]],combine=FALSE)[[1]]
  if(i != 'all'){
    plotlist[[i]] <- plotlist[[i]] + scale_x_reverse(limits=c(21000,0),expand=c(0,0),breaks=seq(0,21000,3000),labels=NULL,minor_breaks=seq(0,21000,1000))+
      labs(x=NULL,y=i)
  }else{
    plotlist[[i]] <- plotlist[[i]] +
      scale_x_reverse(limits=c(21000,0),expand=c(0,0),breaks=seq(0,21000,3000),minor_breaks=seq(0,21000,1000))+
      labs(y=i)
  }
}

ggpubr::annotate_figure(egg::ggarrange(plots=plotlist,ncol=1,padding=0.2),
                top = "Temperature Composite (degC)")



ensOut_all <- list()
for (reg in c(regions$name,'all')){
  #do composite
  if (reg == 'all'){
    fTS_sel <- fTS
  }else{
    fTS_sel <- fTS[pullTsVariable(fTS,'geo_region')==reg]
  }
  ensOut <- compositeEnsembles2(
    fTS = fTS_sel,
    nens = 100,
    binvec = binvec,
    binFun = sampleEnsembleThenBinTs,
    ageVar = "ageEnsemble",
    uncVar = "paleoData_temperature12kUncertainty",
    weights = "geo_weight",
    stanFun = standardizeOverRandomInterval,
    searchRange = c(1000,11000),
    duration = 4000,
    normalizeVariance = TRUE,
    gaussianizeInput = TRUE,
    scale = FALSE,
    spread = TRUE,
    spreadMax=NA,
    samplePct = 1,
  )
  print(reg)
  plot(ensOut)
  ensOut_all[[reg]] <- ensOut
}

lm_model <- lm(y ~ x)

# Extract slope coefficient
cat("Slope of the data:", slope, "\n")



bs <- 1000
slope <- c()
ages <- c()
nens  <- ncol(ensOut$composite)
slopeMatrix <- matrix(NA,nrow=length(ensOut$ages),ncol=nens)
for (age in (ensOut$ages)){
  agei = which(ensOut$ages==age)
  agesSelect <- between(ensOut$ages,age-bs,age+bs)
  if(sum(agesSelect) < bs/100){next}
  #x = as.vector(ensOut$composite[agesSelect,1])[[1]]
  x = apply(ensOut$composite,1,median)[agesSelect]
  y = ensOut$ages[agesSelect]
  ages <- c(ages,age)
  s=coef( lm(x ~ y))["y"]
  print(paste(age,'-',s))
  slope <- c(slope,s*-1000)
  for (i in 1:nens){
    x = as.vector(ensOut$composite[agesSelect,i])[[1]]
    s=coef( lm(x ~ y))["y"]
    slopeMatrix[agei,i] <- s*-1000
  }
}
q = 0.1
ggplot()+
  geom_ribbon(aes(x=ensOut$ages,ymin=apply(slopeMatrix,1,quantile,q/2,na.rm=T),ymax=apply(slopeMatrix,1,quantile,1-q/2,na.rm=T)),
              fill='navy',alpha=0.5)+
  geom_vline(xintercept=c(14700,11700,10500),linetype=2)+
  geom_hline(yintercept=c(0),linetype=1)+
  geom_line(aes(x=ages,y=slope),color='navy')+#,color=slope))+
  scale_x_reverse(limits=c(21000,0),expand=c(0,0),breaks=seq(0,21000,3000),labels=seq(0,21,3),minor_breaks=seq(0,21000,1000))+
  #scale_colour_distiller(palette='RdBu',limits=c(-0.4,0.4),oob=scales::squish)+
  theme_bw()+
  labs(x='age (ka)',y='slope (degC / 1000 yrs)',title='Rate of change in Arctic composite',subtitle=paste0('slope over ',bs*2,' year moving window (shading indicates ',(1-q)*100,'% CI from composite ensemble)'))
  #geom_line(x,slope)+

for (i in 1:ncol(ensOut[["proxyUsed"]])){
  name <- names(ensOut[["proxyUsed"]])[1]
}


Kendall(as.numeric(ensOut$composite[,1]),ensOut$ages)


library(Kendall)

# Generate sample matrix data (replace this with your actual matrix)
set.seed(123)  # for reproducibility
matrix_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # 10x10 matrix

# Function to perform Mann-Kendall test on each row or column of the matrix
perform_mann_kendall <- function(data, dimension = "row") {
  if (dimension == "row") {
    # Apply Mann-Kendall test to each row
    result <- apply(data, 1, function(x) {
      Kendall(x)$p.value
    })
  } else if (dimension == "column") {
    # Apply Mann-Kendall test to each column
    result <- apply(data, 2, function(x) {
      Kendall(x)$p.value
    })
  } else {
    stop("Dimension must be either 'row' or 'column'.")
  }

  return(result)
}

# Perform Mann-Kendall test on rows (assuming rows represent time series)
row_p_values <- perform_mann_kendall(matrix_data, dimension = "row")

# Perform Mann-Kendall test on columns (assuming columns represent time series)
col_p_values <- perform_mann_kendall(matrix_data, dimension = "column")

# Print results
cat("Mann-Kendall test results for rows (time series):\n")
print(row_p_values)
cat("\n")

cat("Mann-Kendall test results for columns (time series):\n")
print(col_p_values)




plotGrid <- function(fTS,color='warmest',col.pal='virids',ageBins=seq(0,12000,400)){
  df <- data.frame(
    tsid = pullTsVariable(fTS,'paleoData_TSid'),
    lat = pullTsVariable(fTS,'geo_latitude'),
    lon = pullTsVariable(fTS,'geo_longitude'),
    cell = pullTsVariable(fTS,'geo_cell'),
    region = pullTsVariable(fTS,'geo_region'),
    weight = pullTsVariable(fTS,'geo_weight'),
    color = NA) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326)
}
  cell <- 26

bs<-1000
binvec<-(seq(0,11000,bs))
binAges<-(seq(0+bs/2,11000,bs))
for (i in 1:nrow(df)){
  ts <- fTS[[which(pullTsVariable(fTS,'paleoData_TSid')==df$tsid[i])]]
  vals <- (simpleBinTs(ts,binvec=binvec))
  df$warmest[i] <- binvec[which.max(vals)]+bs/2
  df$coldest[i] <- binvec[which.min(vals)]+bs/2
}
cellComposites <- list()
count <- 1
for (c in unique(df$cell)[36:45]){
  fTS_sel <- fTS[(pullTsVariable(fTS,'geo_cell')==c)]
  print(paste('cell:',c,'-----',length(fTS_sel),'record(s) -----',
              round(count/length(unique(df$cell))*100),'% of cells'))
  count <- count+1
  if(min(unlist(pullTsVariable(fTS_sel,'age')),na.rm=T)>4000){next}
  ensOut <- compositeEnsembles2(
    fTS = fTS_sel,
    nens = 50,
    binvec = seq(0,11000,200),
    binFun = sampleEnsembleThenBinTs,
    ageVar = "ageEnsemble",
    uncVar = "paleoData_temperature12kUncertainty",
    weights = "geo_weight",
    stanFun = standardizeOverRandomInterval,
    searchRange = c(0,11000),
    duration = 5000,
    normalizeVariance = TRUE,
    gaussianizeInput = FALSE,
    scale = FALSE,
    spread = TRUE,
    samplePct = 1,
  )
  if(sum(!is.na(apply(ensOut$composite,1,median,na.rm=T)))>0){plot(ensOut)}
  cellComposites[[as.character(c)]] <- ensOut
}

polygons$warmest <- NA
polygons$coldest <- NA
for (c in unique(df$cell)){
  ensOut<- cellComposites[[as.character(c)]]
  if(is.null(ensOut)){next}
  composite <- data.frame(ages=binAges,vals=NA)
  for (age in composite$ages){
      idx <- between(ensOut$ages,age-bs/2,age+bs/2)
      composite$vals[composite$ages==age] <- median(apply(ensOut$composite[idx,],2,median,na.rm=T),na.rm=T)
  }
  if(sum(!is.na( composite$vals))>0){
    polygons[polygons$seqnum==c,]$warmest<- composite$age[which.max(composite$vals)]
    polygons[polygons$seqnum==c,]$coldest<- composite$age[which.min(composite$vals)]
  }
}




#Plot polygons
#if (plotcells){
#  cells <- pval2factor(cells,"pvalNet","pvalBinned",color.breaks)
#  map<-map+geom_sf(data=cells, aes(fill=pvalBinned), color=NA)

arcticCircle <- data.frame(id="A",lon=seq(-180,180), lat= rep(58,length(seq(-180,180))))%>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise((geometry = sf::st_combine(geometry))) %>%
  st_cast("MULTILINESTRING")
#get arctic polygons
arctic <- defineRegion(58,90,-180,180)
polygonsArctic  <- suppressWarnings(sf::st_intersection(polygons,arctic)  %>%
                                      mutate(lon=st_coordinates(st_centroid(.))[,1]) %>%
                                      mutate(lat=st_coordinates(st_centroid(.))[,2]) %>%
                                      mutate(area_km2 = units::drop_units(st_area(.)/1000000))) #km^2
#countries to plot
countries <- st_as_sf(rworldmap::getMap("high"))
countriesArctic <- suppressWarnings(sf::st_intersection(countries,arctic))
#Create sf tables
#st_geometry(regions) <- "geometry"
#regions$name <- c("North America","Greenland","Scandinavia","Svalbard","Central Russia","Eastern Russia")
col <- 'warmest'
vals <- polygonsArctic[[col]]
vals <- round(vals/2000)*2000
vals <- as.factor(vals)
color.pal='Purples'
colorvec <- RColorBrewer::brewer.pal(11, 'Spectral')
polygonsArctic$fill <- vals
#map <-
df$warmest<-factor(x=df$warmest,levels=binAges)
polygonsArctic$warmest<-factor(x=polygonsArctic$warmest,levels=binAges)

ggplot() +
  #
    geom_sf(data=polygonsArctic, aes(fill=warmest), color=NA)+
    geom_sf(data=countriesArctic, fill=NA, color='grey30', linewidth=0.07)+
    geom_sf(data=df,aes(fill=warmest),shape=21,size=3)+
    geom_sf(data=regions, fill=NA, color='black', linewidth=0.7)+
  #
    geom_sf(data=arcticCircle, fill=NA, color='black', linewidth=0.7) +
    scale_fill_manual(values=colorvec, drop = FALSE,name='net\np-value',na.value = "white")+
    coord_sf(crs =  CRS("ESRI:102016"),ylim=c(-3300000,3300000),xlim=c(-3300000,3300000)) +
    labs(title=paste("When was the warmest millenium of the Holocene"), alpha='distance\nweight')+
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),        ## <- this line
      panel.border = element_rect(color=NA,fill=NA),
      panel.grid.major = element_blank(), #(color='black',linewidth=0.1),
      # panel.ontop = F
    )

map <- map +
  geom_sf(data=points, aes(fill=pvalBinned),shape=21,size=4) +
  scale_fill_manual(values=colorvec, drop = FALSE,name='net\np-value',na.value = "white")+


  #color.breaks <- c(0.001,0.01,0.025,0.05,0.95,0.975,0.99,0.999)
  colorvec <- RColorBrewer::brewer.pal(length(color.breaks)+1, color.pal)
  colorvec[which(colorvec=="#F7F7F7")]<-'grey'
  if (color.pal=='BrBG'){colorvec<-rev(colorvec)}
  # Point data
  points <- df_raw %>%  filter(time_mid==age) %>%
    st_as_sf(coords = c("geo_longitude", "geo_latitude"), crs = 4326) %>%
    arrange(desc(pvalue_either))
  # Convert pval and bin according to colorvec
  points <- points[unlist(lapply(points$null_probability,sum)) > 0,]
  points$pvalNet <- purrr::map2_dbl(points$pvalue_positive, points$pvalue_negative, handleNet)
  points <- pval2factor(points,"pvalNet","pvalBinned",color.breaks)
  #
  cells <- data.frame()
  plotcells <- TRUE
  if (is.na(poly)){
    plotcells<-FALSE
  } else if (poly == 'cell'){
    cells <- data.frame(cell = sort(unique(points$cell)))
    for (i in 1:nrow(cells)){
      df_cell <- points %>% filter(cell==cells$cell[i])
      if (!("weight" %in% names(df_cell))){df_cell$weight<-NA}
      out <- calculateMultiTestSignificance(df_cell,weights=df_cell$weight,n.ens=median(df_raw$null.hypothesis.n))
      #cells$pvalNet[i] <- out$pvalNet
      cells$pvalNet[i] <- handleNet(out$pvalPos,out$pvalNeg)
    }
    #
    cells <- merge(polygonsArctic,cells,by.x="seqnum",by.y="cell")
  } else if (poly == 'region'){
    cells <- regions
    for (i in 1:nrow(regions)){
      df_cell <- points %>% filter(region==regions$name[i])
      if (nrow(df_cell)==0){
        cells$pvalNet[i]<-NA
        next
      }
      if (!("weight" %in% names(df_cell))){df_cell$weight<-NA}
      out <- calculateMultiTestSignificance(df_cell,weights=df_cell$weight)
      #cells$pvalNet[i] <- out$pvalNet
      cells$pvalNet[i] <- handleNet2(out$pvalPos,out$pvalNeg)
    }
  } else(stop("poly must be 'cell' or 'region'"))
  #
  # Plot
  print((plotcells))
  map <- ggplot() + actR_ggtheme()
    #Plot polygons
  if (plotcells){
    cells <- pval2factor(cells,"pvalNet","pvalBinned",color.breaks)
    map<-map+geom_sf(data=cells, aes(fill=pvalBinned), color=NA)
  }
  map <- map +geom_sf(data=countriesArctic, fill=NA, color='grey30', linewidth=0.07)
  if (!is.null(ice)){
    map <- map + geom_sf(data=ice, fill='white', color='black',linewidth=0.3,alpha=0.3)
  }
  if (plotcells){
    map <- map + geom_sf(data=regions, fill=NA, color='black', linewidth=0.7)
  }
  map <- map +
      geom_sf(data=arcticCircle, fill=NA, color='black', linewidth=0.7) +
      geom_sf(data=points, aes(fill=pvalBinned),shape=21,size=4) +
      scale_fill_manual(values=colorvec, drop = FALSE,name='net\np-value',na.value = "white")+
      coord_sf(crs =  CRS("ESRI:102016"),ylim=c(-3300000,3300000),xlim=c(-3300000,3300000)) +
      labs(title=paste(age/1000, "ka"), alpha='distance\nweight')+
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),        ## <- this line
        panel.border = element_rect(color=NA,fill=NA),
        panel.grid.major = element_blank(), #(color='black',linewidth=0.1),
        # panel.ontop = F
      )


#cell 261, 170,188
#
#
#
#
#
# tsids <- unique(shifts$paleoData_TSid)
#
# df <- data.frame()
# pb <- progress_bar$new(total = length(tsids), format = "[:bar] :percent | ETA: :eta")
# for (i in 1:length(tsids)){
#   sel <- (shifts %>% filter(paleoData_TSid==tsids[i]))[1,28:length(names(shifts))]
#   binmat <- matrix(NA,(length(binyrs)),ncol(sel$time[[1]]))
#   for (c in 1:ncol(binmat)){
#     binmat[,c] <- geoChronR::bin(time = sel$time[[1]][,c],values = sel$paleoData_values[[1]][,c],bin.vec = binvec)[,2]
#   }
#   if (sel$interpDir < 0){
#     ages <- binyrs[unlist(apply(binmat,2,which.min))]
#   }else{
#     ages <- binyrs[unlist(apply(binmat,2,which.max))]
#   }
#   mode <- function(x){
#     ux <- unique(x)
#     m <- ux[which.max(tabulate(match(x, ux)))]
#     return(m)
#   }
#   sel$agesMax <- list(ages)
#   sel$agesMaxMode <- mode(ages)
#   df <- rbind(df,sel)
#   pb$tick()
#
# }
#
# plotMaxValAge <- function(df_raw,age.range,poly=NA, color.breaks = c(0.01,0.05,0.1,0.2,0.8,0.9,0.95,0.99),color.pal='RdBu'){
#   #color.breaks <- c(0.001,0.01,0.025,0.05,0.95,0.975,0.99,0.999)
#   #colorvec <- RColorBrewer::brewer.pal(length(age.range)+2, 'oranges')[2:(length(age.range)+2)]
#   #colorvec[which(colorvec=="#F7F7F7")]<-'grey'
#   # Point data
#   points <- df_raw %>%  st_as_sf(coords = c("geo_longitude", "geo_latitude"), crs = 4326) #%>% arrange(desc(pvalue_either))
#   # Convert pval and bin according to colorvec
#   #points <- points[unlist(lapply(points$null_probability,sum)) > 0,]
#   #points$pvalNet <- purrr::map2_dbl(points$pvalue_positive, points$pvalue_negative, handleNet)
#   #points <- pval2factor(points,"pvalNet","pvalBinned",color.breaks)
#   #
#   cells <- data.frame()
#   plotcells <- TRUE
#   if (is.na(poly)){
#     plotcells<-FALSE
#   } else if (poly == 'cell'){
#     cells <- data.frame(cell = sort(unique(points$cell)),ageSig=NA)
#     for (i in 1:nrow(cells)){
#       df_cell <- points %>% filter(cell==cells$cell[i])
#       if (!("weight" %in% names(df_cell))){df_cell$weight<-NA}
#       ages  <- sort(unique(df_cell$time_mid))
#       pvals <- c() #
#       counts <- c()
#       for (age in ages){
#         df_cell_age <- df_cell %>% filter(time_mid==age)
#         out <- calculateMultiTestSignificance(df_cell_age,weights=df_cell_age$weight,n.ens=median(df_raw$null.hypothesis.n))
#         pvals<- c(pvals, out$pvalEither) #handleNet(out$pvalPos,out$pvalNeg))
#         counts <- c(counts, out$allEventPos)
#       }
#       #If a pvalie tie, go with the one with more positive shifts detected regardless of null
#       if (min(pvals)<0.2){
#         cells$ageSig[i] <- ages[pvals==min(pvals)][which.max(counts[pvals==min(pvals)])]
#       }
#       #cells$pvalNet[i] <- out$pvalNet
#     }
#     #
#     cells <- merge(polygonsArctic,cells,by.x="seqnum",by.y="cell")
#   } else(stop("poly must be 'cell' or 'region'"))
#   #
#   # Plot
#   map <- ggplot() #+ actR_ggtheme()
#   #Plot polygons
#   if (plotcells){
#     #cells$ageSig2 <- as.factor(cells$ageSig)
#
#     #cells <- pval2factor(cells,"pvalNet","pvalBinned",color.breaks)
#     map<-map+geom_sf(data=cells, aes(fill=ageSig2), color=NA)
#   }
#   map <- map + geom_sf(data=countriesArctic, fill=NA, color='grey30', linewidth=0.07)
#   if (plotcells){
#     map <- map + geom_sf(data=regions, fill=NA, color='black', linewidth=0.7)
#   }
#   map<-map +
#     geom_sf(data=arcticCircle, fill=NA, color='black', linewidth=0.7) +
#     geom_sf(data=points, aes(fill=as.factor(agesMaxMode)), color='black',shape=21,size=3) +
#     #scale_fill_manual(values=colorvec, drop = FALSE,name='net\np-value',na.value = "white")+
#     scale_fill_viridis_d(direction=-1)+
#     coord_sf(crs =  sp::CRS("ESRI:102016"),ylim=c(-3300000,3300000),xlim=c(-3300000,3300000)) +
#     #labs(title=paste0(min(age.range/1000),'-',max(age.range/1000),' ka \n age of most significant shift: '))+
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       axis.text = element_blank(),
#       axis.ticks = element_blank(),        ## <- this line
#       panel.border = element_rect(color=NA,fill=NA),
#       panel.grid.major = element_blank(), #(color='black',linewidth=0.1),
#       # panel.ontop = F
#     )
#   return(map)
# }
# plotMaxValAge(df)
#
# z <- (df %>% filter(agesMaxMode<1000) )
# i <- 3
# hist(z[i,]$agesMax[[1]])
# plot(apply(z[i,]$time[[1]],1,median),apply(z[i,]$paleoData_values[[1]],1,median))
