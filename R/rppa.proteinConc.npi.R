rppa.proteinConc.npi <- function(slide, set.neg.to.zero=F, pos.ctrl=c("TOX"), neg.ctrl=c("MIC#1", "MIC#2")){
  require(foreach)
  newSlide <- foreach(currentA = unique(slide$A), .combine = rbind) %dopar% 
  {
    foreach(currentB = unique(slide$B), .combine = rbind) %dopar% 
    {
        foreach(currentFill = unique(slide$Fill), .combine = rbind) %dopar%
        {
          x <- subset(slide, A == currentA & B == currentB & Fill == currentFill)
          y <- x
          y$concentrations <- ((mean(x[x$Sample%in%pos.ctrl, "concentrations"], na.rm=T) - x$concentrations) / 
            (mean(x[x$Sample%in%pos.ctrl, "concentrations"], na.rm=T)-
              mean(x[x$Sample%in%neg.ctrl, "concentrations"], na.rm=T)))  
          y$upper <- (x$upper / x$concentrations) + y$concentrations -1
          y$lower <- (x$lower / x$concentrations) + y$concentrations -1 
          
          if(set.neg.to.zero){
            y[!is.na(y$concentrations) & y$concentrations < 0, "upper"] <- 0
            y[!is.na(y$concentrations) & y$concentrations < 0, "concentrations"] <- 0
          }
          
          return(y)
        }
    }
  }
  
  attr(newSlide, "title") <- attr(slide, "title")
  attr(newSlide, "antibody") <- attr(slide, "antibody")
  attr(newSlide, "blocksPerRow") <- attr(slide, "blocksPerRow")
  
  return(newSlide)
}
