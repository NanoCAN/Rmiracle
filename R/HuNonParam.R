rppa.nonparam <- function(spots, nrep=1, ...){
  require(cobs)
  
  spots <- subset(spots, SpotClass=="Sample")
  
  #convert input table so that each dilution is in one column
  spots.c <- rppa.serialDilution.format(spots, useDepositionsInDilutionSeries=T)
  
  #extract number of different dilutions that are not NA
  numOfDilutions <- attr(spots.c, "numOfDilutions")
  
  #calculate matrix of dilutions
  spots.m <- rppa.serialDilution.dilutionMatrix(spots.c, numOfDilutions, highestDilutionFirst=F)
  
  #noise correction
  #spots.m <- rppa.noiseCorrection(spots.m,highestDilutionFirst=F, subtractNoiseFromAll=T)
  
  nonpa <- getnonpest(spots.m,3)
  
  #combine estimates with signal information
  spots.result <- cbind(spots.c[,1:(ncol(spots.c)-numOfDilutions)], x.weighted.mean=nonpa$x0new, x.err=NA)
  
  spots.summarize <- rppa.serialDilution.summarize(spots.result, useDeposition=F,...)
  spots.summarize$concentrations <- 2^spots.summarize$x.weighted.mean
  spots.summarize$upper <- 0
  spots.summarize$lower <- 0
  
  spots.summarize <- spots.summarize[,!(colnames(spots.summarize) %in% c("x.weighted.mean", "x.err"))]
  attr(spots.summarize, "title") <- attr(spots, "title")
  attr(spots.summarize, "antibody") <- attr(spots, "antibody")
  
  return(spots.summarize)
}