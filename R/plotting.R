

#' Plot an compoite ensemble
#'
#' @param x The output of compositeEnsembles2()
#' @param ... additional parameters (see plotComposite)
#'
#' @return a ggplot2 object
#' @export
plot.paleoComposite <- function(x,...){
  return(plotComposite(x,...))
}






#
#' A function to plot the output of compositeEnsembles2
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom stats median
#' @importFrom geoChronR plotTimeseriesEnsRibbons
#' @importFrom egg ggarrange
#'
#' @param compositeList The output of compositeEnsembles2()
#' @param ageUnits units for the age variable
#' @param valUnits units for the composite (default = standardized)
#' @param title a title for the combined plot (not used if combine == FALSE)
#' @param combine output a combined figure using egg::ggarrange or provide a list of subplots
#' @inheritDotParams geoChronR::plotTimeseriesEnsRibbons probs color.low color.high color.line color.vector line.width
#'
#' @return a combined or list of gg plots
#' @export


plotComposite <- function(compositeList, ageUnits = 'yr BP',valUnits='standardized',title='composite proxy record', combine=TRUE,...){
  ar <- rev(range(compositeList$ages))
  #Plot
  compositeEns <- compositeList$composite[,(apply(!is.na(compositeList$composite),2,sum)>2)]
  if(NCOL(compositeEns)==0){stop('Too many NA. Check results')}
  ensRibbon <- geoChronR::plotTimeseriesEnsRibbons(X=list(values = compositeList$ages, units=ageUnits, variableName='age'),
                                                     Y=list(values = compositeEns,units=valUnits, variableName='composite anomaly'),...)+
      theme_bw()
  # Add x axis
  ensRibbon <- suppressMessages(ensRibbon +
                                  scale_x_reverse(limits=ar, expand=c(0,0),
                                                  name=ifelse(combine, ' ', paste0('age (',ageUnits,')')),
                                                  position=ifelse(combine, 'top', 'bottom'))
  )

  #
  #Plot Count
  bs = abs(stats::median(diff(compositeList$ages)))
  tsAvailability <- ggplot2::ggplot()+
    geom_bar(stat="identity",aes(x=compositeList$ages,y=apply(compositeList$proxyUsed>0,1,sum)),
             fill='lightgrey',color='darkgrey',width=bs)+
    scale_y_continuous(limits=c(0,ncol(compositeList$proxyUsed)*1.1),expand=c(0,0),name='count')+
    scale_x_reverse(limits=c(ar[1]+bs/2,ar[2]-bs/2),expand=c(0,0),name=paste0('age (',ageUnits,')')) + theme_bw()
  #
  #Combine
  plots <- list(compositeEns=ensRibbon,dataUsed=tsAvailability)
  if (combine){
    plots <-   egg::ggarrange(plots=plots,ncol=1,heights=c(2,1),top=title)
  }
  return(plots)
}







#' Print composite output
#'
#' @param x composite output
#' @param ... additional inputs (see printShift)
#' @export
print.paleoComposite <- function(x,...){
  printComposite(x,...)
}


#' Print compositeEnsembles2() output
#' @importFrom stats median
#' @importFrom glue glue
#' @import dplyr
#'
#' @param compositeList The output of compositeEnsembles2()
#' @export
printComposite <- function(compositeList){ #com params.to.print = c("cpt.fun","minimum.segment.length","method","penalty","ncpts.max")){
  #
  new_df <- data.frame(age=compositeList$ages,
                       median=round(apply(compositeList$composite,1,median,na.rm=T),4),
                       count=round(apply(compositeList$proxyUsed,1,sum,na.rm=T)))
  messages <- list()

  #Composote created
  ar <- crayon::green(paste(range(new_df$age), collapse='-'))
  res <-  crayon::green(median(diff(new_df$age)))
  messages[[1]] <- glue("Composite created for ages: {ar} yr BP at {res} year resolution")

  #Number of records
  proxyN <- crayon::green(NCOL(compositeList$proxyUsed))
  ar <- crayon::green(paste(range((new_df %>% dplyr::filter(count==max(new_df$count,na.rm=T)))$age), collapse='-'))
  messages[[2]] <- glue("{proxyN} proxy records used. The maximum data density occurs between {ar} yr BP")

  #ensembles
  messages[[3]] <- glue("{crayon::green(NCOL(compositeList$composite))} ensembles generated")
  agemin <- new_df$age[new_df$median==min(new_df$median)]
  agemax <- new_df$age[new_df$median==max(new_df$median)]
  messages[[4]] <- glue("The maximum median value occurs at {crayon::red(agemax)} yr BP.")
  messages[[5]] <- glue("The minimum median value occurs at {crayon::blue(agemin)} yr BP.")

  #print
  for (m in messages){print(m)}
  print(utils::head(new_df%>%filter(!is.na(median))))
  #print(new_df)
}
