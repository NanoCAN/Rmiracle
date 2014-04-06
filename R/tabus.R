rppa.tabus <- function(spots, nrep=1, ...){
  
  spots$Signal <- spots$FG-spots$BG #method can't deal with NA values
  spots$Signal[spots$Signal < 0] <- 0
  
  spots <- subset(spots, SpotClass=="Sample")
  
  #convert input table so that each dilution is in one column
  spots.c <- rppa.serialDilution.format(spots)
  
  #extract number of different dilutions that are not NA
  numOfDilutions <- length(unique(spots$DilutionFactor[!is.na(spots$DilutionFactor)]))
  #calculate matrix of dilutions
  spots.m <- rppa.serialDilution.dilutionMatrix(spots.c, numOfDilutions, highestDilutionFirst=F)
  tabus<-dilutionFitrep(spots.m, nrep);
  #fittedData <- data.frame(x=spots$Signal, y=predict(tabus$fit, data.frame(x=spots$Signal)))
  
  #combine estimates with signal information
  spots.result <- cbind(spots.c[,1:(ncol(spots.c)-numOfDilutions)], x.weighted.mean=tabus$pass2, x.err=NA)
  
  spots.summarize <- rppa.serialDilution.summarize(spots.result, ...)
  spots.summarize$concentrations <- 2^spots.summarize$x.weighted.mean
  spots.summarize$upper <- spots.summarize$concentrations
  spots.summarize$lower <- spots.summarize$concentrations
  
  spots.summarize <- spots.summarize[,!(colnames(spots.summarize) %in% c("x.weighted.mean", "x.err"))]
  attr(spots.summarize, "title") <- attr(spots, "title")
  attr(spots.summarize, "antibody") <- attr(spots, "antibody")
  
  readout <- data.frame(concentrations=spots.summarize$readout,
                        upper=spots.summarize$readout.sem + spots.summarize$readout,
                        lower=spots.summarize$readout - spots.summarize$readout.sem)
  readout.centered <- readout / mean(readout$concentrations, na.rm=T)
  attr(spots.summarize, "readout") <- readout
  attr(spots.summarize, "readout.centered") <- readout.centered
    
  return(spots.summarize)
}