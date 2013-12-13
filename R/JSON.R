rppa.load <- function (connection=NULL, barcode=NA, slideIndex=NA, securityToken=NA, baseUrl = "http://localhost:8080/MIRACLE/spotExport/", filter.bad.signals=T, apply.shifts=T) 
{
  require(RJSONIO)
  require(plyr)
  
  #no means of authentication given
  if(is.null(connection) && is.na(securityToken)){
    cat("Cannot authenticate. Either login using rppa.authenticate or provide a security token\n")
  }
  
  #authentication via securityToken
  if(is.null(connection) && !is.na(securityToken)){
    cat("Using security token authentication\n")
    
    #create curl handle without authentication
    agent="Mozilla/5.0" #or whatever 
    
    #Set RCurl pars
    connection = getCurlHandle()
    curlSetOpt(ssl.verifypeer=FALSE, timeout=60, cookiefile=tempfile(), cookiejar=tempfile(), useragent = agent, followlocation = TRUE, curl=connection, verbose=FALSE)
    
    #first check if security token is valid
    isTokenValid <- scan(text=getURL(paste(baseUrl, "isSecurityTokenValid/", securityToken, sep = ""), curl=connection), what="boolean")
    if(isTokenValid) slideIndex <- scan(text=getURL(paste(baseUrl, "getSlideIdFromSecurityToken/", securityToken, sep = ""), curl=connection), what="integer")
    else{
      cat("Security token is invalid!\n")
      return(NA)
    } 
  }
  
  else{
    cat("Using user authentication\n")
  }
  
  #check input
  if(is.na(slideIndex) && is.na(barcode)){
    cat("You have to specify either barcode or slideIndex")
    return(NA)
  } 
  
  else if(is.na(slideIndex))
  {
    barcodeUrl <- paste(baseUrl, "getIdFromBarcode/", barcode, sep = "")
    if(!is.na(securityToken)) barcodeUrl <- paste(barcodeUrl, "?securityToken=", securityToken, sep="")
    slideIndex <- getURL(barcodeUrl, curl=connection)
    slideIndex <- fromJSON(slideIndex)
    if(length(slideIndex) > 1)
    {
      cat(paste("There are several options for barcode", barcode,". Please make a selection:"))
      for(i in slideIndex){
        titleUrl <- paste(baseUrl, "getTitle/", i, sep = "")
        if(!is.na(securityToken)) titleUrl <- paste(titleUrl, "?securityToken=", securityToken, sep="")
                           
        cat(paste(i, ":", paste(scan(text=getURL(titleUrl, curl=connection), what="character"), collapse=" ")))  
      }
      slideIndex <- readline("Please enter a slide index:")
    }
    cat(paste("Barcode", barcode, "was identified as slide index", slideIndex, "\n"), sep=" ")
  }
  
  #read the data from database
  cat(paste("Reading spots for slide index", slideIndex, "\n"))
  spotsUrl <- paste(baseUrl, "exportAsJSON/", slideIndex, sep = "")
  if(!is.na(securityToken)) spotsUrl <- paste(spotsUrl, "?securityToken=", securityToken, sep="")
  
  spots <- getURL(spotsUrl, curl=connection)
  spots <- ldply(fromJSON(spots, simplify = T, nullValue = NA))
  cat(paste(dim(spots)[1], "spots read. Formatting...\n"))
  
  metaUrl <- paste(baseUrl, "exportMetaDataAsJSON/", slideIndex, sep = "")
  if(!is.na(securityToken)) metaUrl <- paste(metaUrl, "?securityToken=", securityToken, sep="")
  
  meta <- getURL(metaUrl, curl=connection)
  meta <- fromJSON(meta, simplify=T)
  
  colnames(spots) <- meta
  
  #reformat column types
  spots$id <- as.integer(spots$id)
  spots$FG <- as.double(spots$FG)
  spots$BG <- as.double(spots$BG)
  spots$Signal <- as.double(spots$Signal)
  spots$Diameter <- as.double(spots$Diameter)
  spots$Flag <- as.double(spots$Flag)
  spots$Block <- as.integer(spots$Block)
  spots$Row <- as.integer(spots$Row)
  spots$Column <- as.integer(spots$Column)
  spots$CellLine <- as.factor(spots$CellLine)
  spots$Treatment <- as.factor(spots$Treatment)
  spots$Inducer <- as.factor(spots$Inducer)
  spots$LysisBuffer <- as.factor(spots$LysisBuffer)
  spots$SampleName <- as.factor(spots$SampleName)
  spots$SampleType <- as.factor(spots$SampleType)
  spots$TargetGene <- as.factor(spots$TargetGene)
  spots$DilutionFactor <- as.double(spots$DilutionFactor)
  spots$PlateCol <- as.integer(spots$PlateCol)
  spots$PlateRow <- as.integer(spots$PlateRow)
  spots$PlateLayout <- as.integer(spots$PlateLayout)
  
  #add shifts
  shiftUrl <- paste(baseUrl, "exportShiftsAsJSON/", slideIndex, sep = "")
  if(!is.na(securityToken)) shiftUrl <- paste(shiftUrl, "?securityToken=", securityToken, sep="")
  
  shifts <- getURL(shiftUrl, curl=connection)
  if (length(fromJSON(shifts)) > 0) {
    shifts <- ldply(fromJSON(shifts, simplify = T, nullValue = NA))
    spots <- merge(spots, shifts, by = "Block", all.x = T)
  }
  else{
    spots$vshift <- 0
    spots$hshift <- 0
  }
  #apply shifts
  spots <- rppa.vshift(spots)
  spots <- rppa.hshift(spots)
  
  #add depositions
  deposUrl <- paste(baseUrl, "getDepositionPattern/", slideIndex, sep = "")
  if(!is.na(securityToken)) deposUrl <- paste(deposUrl, "?securityToken=", securityToken, sep="")
  
  depositionPattern <- scan(text=getURL(deposUrl, curl=connection), what = "integer")
  depositionPattern <- gsub("\\[", "", depositionPattern)
  depositionPattern <- gsub("\\]", "", depositionPattern)
  depositionPattern <- as.integer(strsplit(depositionPattern, 
                                           ",")[[1]])
  spots$Deposition <- spots$Column%%length(depositionPattern)
  spots$Deposition[spots$Deposition == 0] <- length(depositionPattern)
  spots$Deposition <- depositionPattern[spots$Deposition]
  
  #filter bad signals
  if(filter.bad.signals)
  {
    spots <- rppa.filter.diameter(spots)
    spots <- rppa.filter.flagged(spots)
    spots <- rppa.filter.neg.values(spots)
  }
  blocksUrl <- paste(baseUrl, "getBlocksPerRow/", slideIndex, sep = "")
  if(!is.na(securityToken)) blocksUrl <- paste(blocksUrl, "?securityToken=", securityToken, sep="")
  
  titleUrl <- paste(baseUrl, "getTitle/", slideIndex, sep = "")
  if(!is.na(securityToken)) titleUrl <- paste(titleUrl, "?securityToken=", securityToken, sep="")
  
  antibodyUrl <- paste(baseUrl, "getAntibody/", slideIndex, sep = "")                     
  if(!is.na(securityToken)) antibodyUrl <- paste(antibodyUrl, "?securityToken=", securityToken, sep="")
  
  spots <- rppa.set.blocksPerRow(spots, as.integer(scan(text=getURL(blocksUrl, curl=connection), what = "integer")))
  spots <- rppa.set.title(spots, paste(scan(text=getURL(titleUrl, curl=connection), what = "character"), collapse=" "))
  spots <- rppa.set.antibody(spots, paste(scan(text=getURL(antibodyUrl, curl=connection), what = "character"), collapse=" "))
  cat("...everything done. returning data.")
  return(spots)
}

rppa.filter.diameter <- function(spots)
{
  spots$Diameter <- as.double(spots$Diameter)
  spots$Signal[spots$Diameter >= 250] <- NA
  return(spots)
}

rppa.filter.neg.values <- function(spots)
{
  spots$FG <- as.double(spots$FG)
  spots$BG <- as.double(spots$BG)
  spots$Signal[(spots$FG-spots$BG <= 0)] <- NA
  return(spots)
}

rppa.filter.flagged <- function(spots)
{
  spots$Flag <- as.double(spots$Flag)
  spots$Signal[spots$Flag != 0] <- NA
  return(spots)
}