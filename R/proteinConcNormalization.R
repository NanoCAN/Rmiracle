rppa.proteinConc.normalize <- function(slideA, slideB, method="houseKeeping", normalize.with.median.first = T, 
                                       target.column="Slide", normalize.per.deposition=F, output.all=F)
{ 
  #which index is the median value selected for normalization. 
  which.median = function(x) {
    if (length(x) %% 2 != 0) {
      which(x == median(x))
    } else if (length(x) %% 2 == 0) {
      a = sort(x)[c(length(x)/2, length(x)/2+1)]
      c(which(x == a[1]), which(x == a[2]))
    }
  }
  
  #check if supported method has been selected
  if(!(method %in% c("houseKeeping", "medianLoading", "variableSlope"))){
    cat("specify a method. Either houseKeeping or medianLoading/variableSlope if a list of slides is provided")
    return(NA)
  }
  
  if(is.data.frame(slideB)) slideB <- list(slideB)
  numberOfSlidesInB = length(slideB)
  
  #method for median normalization
  sub.normalize <- function(slide){
    slide$upper <- slide$upper / median(slide$concentrations, na.rm=T)
    slide$lower <- slide$lower / median(slide$concentrations, na.rm=T)
    slide$concentrations <- slide$concentrations / median(slide$concentrations, na.rm=T)
    return(slide)
  }
  
  #first median normalization for target slide
  if(normalize.with.median.first)
  {
    if(!normalize.per.deposition){
      slideA <- sub.normalize(slideA)
    }
    else{
      slideA <- ddply(slideA, .(Deposition), sub.normalize)
    }
  }
  
  #save target slides confidence intervals as percentage, so we can calculate them later
  #based on corrected values (variable slope normalization)
  slideA$upper <- (( slideA$upper - slideA$concentrations) /slideA$concentrations )
  slideA$lower <- (( slideA$concentrations - slideA$lower) / slideA$concentrations )
  
  
  #normalize with median
  if(normalize.with.median.first)
  {
    if(!normalize.per.deposition){
      slideB <- lapply(slideB, sub.normalize)
    }
    else{
      slideB <- foreach(slide=slideB) %do% ddply(slide, .(Deposition), sub.normalize)
    }
  }
    
  #Firstly, we collect all necessary values and bind the columns (one column per slide)
  allConcentrations <- foreach(slide=slideB, .combine=cbind) %do%{
    return(slide$concentrations)
  }
    
  #confidence intervals are kept as percentage of their concentration estimates, so
  #we can reuse them after applying a correction factor.
  allUpper <- foreach(slide=slideB, .combine=cbind) %do%{
    return(( slide$upper - slide$concentrations) / slide$concentrations )
  }
  
  allLower <- foreach(slide=slideB, .combine=cbind) %do%{
    return((slide$concentrations - slide$lower) /slide$concentrations )
  }

  slideB <- slideB[[1]]
  if(numberOfSlidesInB > 1)
  {
    #we choose to keep the first one for its layout properties and overwrite its title
    attr(slideB, "antibody") <- "mixed"
    attr(slideB, "title") <- method
  }
 
  if(method == "variableSlope")
  {
    #apply to all slides, including target slide
    allSlideConcentrations <- cbind(slideA$concentrations, allConcentrations)
    gammaB <- estimateGamma(allSlideConcentrations,method="other")
    #exchanged '-' and '/' since we don't use log values
    
    mydata.gam <- sweep(allSlideConcentrations,2,gammaB,"/")
    slideA$concentrations <- mydata.gam[,1]
    allConcentrations <- mydata.gam[,-1]
  }
  if(method %in% c("medianLoading", "variableSlope") && numberOfSlidesInB > 1)
  {
    #in median loading normalization the median value of all housekeeping proteins is selected.
    #variable slope normalization follows the same principle, but applies a correction factor gamma first.
    
    #take the median of each row
    medians <- apply(allConcentrations, 1, which.median)
    
    #index of the median is needed for finding the 
    #corresponding confidence intervals
    for(i in 1:length(medians)){
      slideB$concentrations[i] <- allConcentrations[i, medians[i]]
      slideB$upper[i] <- allUpper[i, medians[i]] 
      slideB$lower[i] <- allLower[i, medians[i]] 
    }
  }
  
  poolError <- function(errorRates){
    apply(errorRates, 1, function(x){
      return(sqrt(sum(sapply(x, function(y){return(y*y)}))) / length(x))
    })
  }
  
  poolDivisionError <- function(errorRates){
    apply(errorRates, 1, function(x){
      return(sqrt(sum(sapply(x, function(y){return(y*y)}))))
    })
  }
  
  if(method == "houseKeeping" && numberOfSlidesInB > 1){
    slideB$concentrations <- apply(allConcentrations, 1, mean)
    slideB$upper <- poolError(allUpper)
    slideB$lower <- poolError(allLower)
  }
  else if(method == "houseKeeping"){
    slideB$upper <- allUpper
    slideB$lower <- allLower
  }
  
  #check if target column is free on both slides
  if(output.all && !is.null(slideA[[target.column]]))
  {
    cat("The target column of slideA is not available.")
    return(NA)
  }  
  if(output.all && !is.null(slideB[[target.column]]))
  {
    cat("The target column of slideB is not available.")
    return(NA)
  }
  slideA[[target.column]] <- attr(slideA, "title")
  slideB[[target.column]] <- attr(slideB, "title")
  
  result <- slideA
  result$concentrations <- slideA$concentrations / slideB$concentrations
  
  #combine error
  result$upper <- poolDivisionError(cbind(slideA$upper,slideB$upper)) 
  result$lower <- poolDivisionError(cbind(slideA$lower,slideB$lower))
  
  result$upper <- (result$upper * result$concentrations) + result$concentrations
  result$lower <- result$concentrations - (result$lower * result$concentrations)
  
  result[[target.column]] <- paste(slideA[[target.column]], "normalized by", slideB[[target.column]])
  if(output.all) result <- rbind(slideA, result, slideB)
    
  return(result)
}

rppa.specific.dilution <- function(spots, dilution=0.25, deposition=4, ...)
{
  #if(!is.null(spots$Inducer)) spots$Inducer <- gsub(" [0-9]+[.][0-9] mM", "", spots$Inducer )
  
  spots.subset <- subset(spots, DilutionFactor == dilution & Deposition == deposition & !is.na(Signal))
  spots.subset$x.weighted.mean <- spots.subset$Signal
  spots.subset$x.err <- 0
  spots.summarize <- rppa.serialDilution.summarize(spots.subset, ...)
  spots.summarize$x.err <- spots.summarize$sem
  
  if(length(spots.summarize$x.err[is.na(spots.summarize$x.err)])>0) cat("WARNING: some samples have only one valid value and no standard error could be computed!")
  
  
  spots.summarize$concentrations <- spots.summarize$x.weighted.mean
  
  spots.summarize$upper <- spots.summarize$x.weighted.mean + spots.summarize$x.err
  spots.summarize$lower <- spots.summarize$x.weighted.mean - spots.summarize$x.err
  
  spots.summarize <- spots.summarize[,!(colnames(spots.summarize) %in% c("sem", "x.weighted.mean", "x.err", "Deposition"))]
  
  attr(spots.summarize, "title") <- attr(spots, "title")
  attr(spots.summarize, "antibody") <- attr(spots, "antibody")
  
  return(spots.summarize)
}

rppa.duplicate.nas <- function(data.protein.conc.copy)
{
  foreach(property=c("A","B")) %do%{
    data.protein.conc.copy <- foreach(i=1:nrow(data.protein.conc.copy), .combine=rbind) %do%  {
      if(is.na(data.protein.conc.copy[i,property])){
        foreach(A=levels(data.protein.conc.copy[[property]]), .combine=rbind) %do% {
          currentRow <- data.protein.conc.copy[i,]
          currentRow[[property]] <- A
          return(currentRow)
        }
      }
      else return(data.protein.conc.copy[i,])
    }
  }
  return(data.protein.conc.copy)
}

rppa.normalize.to.ref.sample <- function(data.protein.conc, sampleReference, each.A=F, each.B=F, specific.A=NULL, specific.B=NULL, each.fill=F, method="mean")
{
  library(plyr)
  
  toRefSample <- function(data.protein.conc){ 
    if(is.null(specific.A) && is.null(specific.B)) my.subset <- subset(data.protein.conc, Sample %in% sampleReference)
    else if(is.null(specific.B)){
      if(!is.null(specific.A)) my.subset <- subset(data.protein.conc, Sample %in% sampleReference & A == specific.A)
      else my.subset <- subset(data.protein.conc, Sample %in% sampleReference & is.na(A))
    } 
    else if(is.null(specific.A)){
      if(!is.null(specific.B)) my.subset <- subset(data.protein.conc, Sample %in% sampleReference & B == specific.B)
      else  my.subset <- subset(data.protein.conc, Sample %in% sampleReference & is.na(B))
    } 
    else{
      if(!is.null(specific.A) && !is.null(specific.B)) my.subset <- subset(data.protein.conc, Sample %in% sampleReference & A == specific.A & B == specific.B)
      else if(is.null(specific.A) & !is.null(specific.B)) my.subset <- subset(data.protein.conc, Sample %in% sampleReference & B == specific.B & is.na(A))
      else if(is.null(specific.B) & !is.null(specific.A)) my.subset <- subset(data.protein.conc, Sample %in% sampleReference & A == specific.A & is.na(B))
      else my.subset <- subset(data.protein.conc, Sample %in% sampleReference & is.na(A) & is.na(B))  
    } 
       
    if(method == "mean")  meanOfRefSample <- mean(my.subset$concentrations, na.rm=T)
    else if(method == "median")  meanOfRefSample <- median(my.subset$concentrations, na.rm=T)
       
    data.protein.conc <- within(data.protein.conc, {
      concentrations <- concentrations / meanOfRefSample  
      upper <- upper / meanOfRefSample
      lower <- lower / meanOfRefSample
    }, meanOfRefSample=meanOfRefSample)
  }
  if(each.fill)
  {
    data.protein.conc <- ddply(data.protein.conc, intersect(colnames(data.protein.conc), c("A", "B", "Slide", "Fill")), function(x, sampleRef){ 
      within(x, {
              reference <- mean(concentrations[Sample %in% sampleRef], na.rm=T)
              concentrations <- concentrations / reference
              upper <- upper / reference
              lower <- lower / reference
           })
    }, sampleRef=sampleReference)
    data.protein.conc <- data.protein.conc[,setdiff(colnames(data.protein.conc), c(".id", "reference"))]
  }
  
  else if(each.A && each.B){  
    data.protein.conc <- ddply(data.protein.conc, .(A, B), toRefSample)
  }
  else if(each.A){
    data.protein.conc <- ddply(data.protein.conc, .(A), toRefSample)
  }
  else if(each.B){
    data.protein.conc <- ddply(data.protein.conc, .(B), toRefSample)
  }
  else 
  {
    data.protein.conc <- toRefSample(data.protein.conc)
  }
  
  return(data.protein.conc)
}