rppa.quantile.normalize <- function(sdcs, method){
  concentrations <- foreach(sdc=sdcs, .combine=cbind) %do% {
    sdc$concentrations
  }
  corrected.values <- concentrations
  
  if(method=="aromaQuant")
  {
    library(aroma.light)  
    corrected.values <- normalizeQuantileSpline(as.matrix(concentrations))    
  }
  else if(method=="affyQuant")
  {
    library(affy)
    corrected.values <- normalize.qspline(concentrations, samples=4, na.rm=T)
  }
  else if(method=="cyclicLoess")
  {
    library(limma)
    corrected.values <- normalizeCyclicLoess(concentrations)
  }
  else{
    stop("you have to select a method")
  }
  for(i in 1:length(sdcs)){
    sdc <- sdcs[[i]]
    sdc$upper <- (( sdc$upper - sdc$concentrations) / sdc$concentrations )
    sdc$lower <- (( sdc$concentrations - sdc$lower) / sdc$concentrations )
    sdc$concentrations <- corrected.values[,i] 
    sdc$upper <- (sdc$upper * sdc$concentrations) + sdc$concentrations
    sdc$lower <- sdc$concentrations - (sdc$lower * sdc$concentrations)
    sdcs[[i]] <- sdc
  }
  return(sdcs)
}