rppa.load.all.pmt <- function(connection, barcode, baseUrl="http://10.149.64.36:8080/MIRACLE/spotExport/")
{
  library(foreach)
  library(RCurl)
  
  barcodeUrl <- paste(baseUrl, "getIdFromBarcode/", barcode, sep = "")
  slideIndex <- getURL(barcodeUrl, curl=connection)
  slideIndex <- fromJSON(slideIndex)
  
  pmts <- foreach(slide=slideIndex, .combine=cbind) %do% {
    pmtUrl <- paste(baseUrl, "getPMT/", slide, sep = "")   
    as.integer(scan(text=getURL(pmtUrl, curl=connection), what = "integer"))
  } 
  slideIndex = slideIndex[order(pmts)]
  
  slides <- foreach(slide=slideIndex) %do% {
    rppa.load(connection=connection,slideIndex=slide,baseUrl=baseUrl)
  }
  return(slides)
}

rppa.pick.best.pmt <- function(slides)
{
  quality <- foreach(slide=slides, .combine=rbind) %do% {
    slide_iqr <- IQR(na.omit(subset(slide, SpotClass="Sample")$FG))
    
    spots <- subset(slide, SpotClass=="Sample")
    
    #convert input table so that each dilution is in one column
    spots.c <- rppa.serialDilution.format(spots)
    
    #extract number of different dilutions that are not NA
    numOfDilutions <- length(unique(spots$DilutionFactor[!is.na(spots$DilutionFactor)])) 
    
    #calculate matrix of dilutions
    spots.m <- rppa.serialDilution.dilutionMatrix(spots.c, numOfDilutions)
    slide_dil_corr <- mean(cor(spots.m,use="pairwise.complete.obs"))
    
    result <- data.frame(iqr=slide_iqr, dilution_correlation=slide_dil_corr)
    return(result)
  }
  print(quality)
  return(slides[[which.max(quality[,1] * quality[,2])]])
}

rppa.multiscan.aroma <- function(slides){
  library(foreach)
  library(aroma.light)
  library(plyr)
  slides <- foreach(slide=slides) %do% subset(slide, SpotClass=="Sample")
  signal <- foreach(slide=slides, .combine=cbind) %do% slide$FG
  fit <- calibrateMultiscan(signal, center=FALSE, satSignal=55000)
  slide <- slides[[1]]
  attr(slide, "PMT") <- "mixed"
  slide$Signal <- as.numeric(fit)
  return(slide)
}

rppa.multiscan.lyng <- function(slides, highPMT, lowPMT){
  highPMT <- slides[[highPMT]]
  lowPMT <- slides[[lowPMT]]
  saturated <- highPMT$FG > 55000 & highPMT$SpotClass == "Sample"
  selected <- highPMT$FG > 20000 & highPMT$FG < 50000
  
  if(sum(saturated, na.rm=T) == 0) stop("no saturated spots found in higher PMT")
  if(sum(lowPMT[lowPMT$SpotClass == "Sample", "FG"] > 50000, na.rm=T) > 0) stop("lower PMT may not contain saturated spots")
  #calculate correction factor
  F_cor <- highPMT[selected, "Signal"] / lowPMT[selected, "Signal"]
  F_cor <- sum(F_cor, na.rm=T) / length(na.omit(F_cor))
  cat(F_cor)
  
  #apply correction factor to saturated spots
  highPMT[saturated, "Signal"] <- lowPMT[saturated, "Signal"] * F_cor
  
  return(highPMT)
}
