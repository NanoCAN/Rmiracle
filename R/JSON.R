rppa.load.readout <- function(connection=NULL, securityToken=NA, readoutIndex=NA, baseUrl = "http://localhost:8080/MIRACLE/readoutExport/"){
  
  library(RJSONIO)
  library(plyr)
  library(Hmisc)
  
  connection <- rppa.check.connection(connection, securityToken, baseUrl)  
  if(!is.na(securityToken) && is.na(readoutIndex)) readoutIndex <- scan(text=getURL(paste(baseUrl, "getReadoutIdFromSecurityToken/", securityToken, sep = ""), curl=connection), what="integer")
  
  #read the data from database
  cat(paste("Fetching readout for id", readoutIndex, "\n"))
  readoutUrl <- paste(baseUrl, "exportAsJSON/", readoutIndex, sep = "")
  if(!is.na(securityToken)) readoutUrl <- paste(readoutUrl, "?securityToken=", securityToken, sep="")
  
  wells <- getURL(readoutUrl, curl=connection)
  wells <- ldply(fromJSON(wells, simplify = T, nullValue = NA))
  cat(paste(dim(wells)[1], "well readouts downloaded. Formatting...\n"))
  
  metaUrl <- paste(baseUrl, "exportMetaDataAsJSON/", readoutIndex, sep = "")
  if(!is.na(securityToken)) metaUrl <- paste(metaUrl, "?securityToken=", securityToken, sep="")
  
  meta <- getURL(metaUrl, curl=connection)
  meta <- fromJSON(meta, simplify=T)
  colnames(wells) <- meta
  wells <- rppa.reformatPlateColTypes(wells)

  assayUrl <- paste(baseUrl, "getAssayType/", readoutIndex, sep = "")
  if(!is.na(securityToken)) assayUrl <- paste(assayUrl, "?securityToken=", securityToken, sep="")
  wells$AssayType <- paste(scan(text=getURL(assayUrl, curl=connection), what = "character"), collapse="")
  
  return(wells)
}

rppa.batch.load.readouts <- function(connection=NULL, readoutSecurityTokens=NULL, plateSecurityTokens=NULL, readoutIndices=NULL, baseUrl = "http://localhost:8080/MIRACLE/readoutExport/"){
  library(foreach)
  if(is.null(readoutSecurityTokens) &&is.null(plateSecurityTokens) && is.null(readoutIndices)) stop("you have to provide securityTokens or indices together with a connection object.")
  
  if(is.null(readoutSecurityTokens) && !is.null(plateSecurityTokens))
  {
    connection <- rppa.check.connection(connection, plateSecurityTokens[1], baseUrl)  
    
    readoutSecurityTokens <- foreach(token=plateSecurityTokens, .combine=cbind) %do% {
      plateTokens <- getURL(paste(baseUrl, "getReadoutSecurityTokensFromPlateSecurityToken/", token, sep = ""), curl=connection)
      plateTokens <- fromJSON(plateTokens, simplify = T, nullValue = NA)
      if(length(plateTokens) > 1) stop("Rmiracle currently supports only one readout per plate in the analysis")
      else return(plateTokens)
    }
  }
  
  if(!is.null(readoutSecurityTokens)){
    result <- foreach(token=readoutSecurityTokens, .combine=rbind) %do%{
      readout <- rppa.load.readout(connection, token, baseUrl=baseUrl)
    } 
  }
  else if(!is.null(connection) && !is.null(readoutIndices)){
    result <- foreach(index=readoutIndices, .combine=rbind) %do% rppa.load.readout(connection, readoutIndex=index, baseUrl=baseUrl)
  }
  
  return(result)
}

rppa.check.connection <- function(connection=NULL, securityToken=NA, baseUrl){
  library(RCurl)
  
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
    isTokenValid <- scan(text=capitalize(getURL(paste(baseUrl, "isSecurityTokenValid/", securityToken, sep = ""), curl=connection)), what=TRUE)
    if(!isTokenValid){
      stop("Security token is invalid!\n")
    } 
  }
  else{
    message("Using user authentication\n")
  }
  return(connection)
}

#' Load an RPPA slide from MIRACLE
#'
#' This function loads a RPPA slide from MIRACLE using one of the 
#' following as input: a connection object obtained with \code{\link{rppa.authenticate}}
#' together with either barcode or slideIndex or a securityToken where user authentication
#' is completly omitted. In addition we need to specify MIRACLE's URL as baseUrl, specify signal
#' filter options and optionally correct for block shifts in the process.
#'
#' @param connection RCurl connection object obtained through \code{\link{rppa.authenticate}}
#' @param barcode Barcode of the slide
#' @param slideIndex Database ID of the slide in MIRACLE
#' @param securityToken UUID of the slide used for secure access without \code{connection}
#' @param baseUrl URL pointing to the spotExport controller of the used MIRACLE instance
#' @param filter.diameter if true filter spots with diameter > 200 (set to NA)
#' @param filter.negative if true filter negative values (set to NA)
#' @param filter.flag if true filter flagged values, e.g. where flag property is not 0 (set to NA)
#' @param apply.shifts if true apply horizontal and vertical shifts
#' @keywords json io load read
#' @export
#' @examples
#' baseUrl <- "http://192.168.0.1:8080/MIRACLE/spotExport/"
#' conn <- rppa.authenticate(user="mlist", password="xxx", baseUrl=baseUrl)
#' rppa.load(connection=conn, slideIndex=1, baseUrl=baseUrl)
#' rppa.load(connection=conn, barcode="ABCD")
#' rppa.load(securityToken="abcd-xzft-djfk-2345")
rppa.load <- function (connection=NULL, barcode=NA, slideIndex=NA, securityToken=NA, baseUrl = "http://localhost:8080/MIRACLE/spotExport/", filter.diameter=T, filter.neg.values=T, filter.flag=T, apply.shifts=T) 
{
  library(RJSONIO)
  library(plyr)
  library(Hmisc)
  
  connection <- rppa.check.connection(connection, securityToken, baseUrl)
  
  if(!is.na(securityToken) && is.na(slideIndex)) slideIndex <- scan(text=getURL(paste(baseUrl, "getSlideIdFromSecurityToken/", securityToken, sep = ""), curl=connection), what="integer")
 
  #check input
  if(is.na(slideIndex) && is.na(barcode)){
    stop("You have to specify either barcode or slideIndex")
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
  
  spots <- rppa.reformatColTypes(spots)

  #add shifts
  shiftUrl <- paste(baseUrl, "exportShiftsAsJSON/", slideIndex, sep = "")
  if(!is.na(securityToken)) shiftUrl <- paste(shiftUrl, "?securityToken=", securityToken, sep="")
  
  shifts <- getURL(shiftUrl, curl=connection)
  if(length(fromJSON(shifts)) > 0) {
    shifts <- ldply(fromJSON(shifts, simplify = T, nullValue = NA))
    spots <- merge(spots, shifts, by = "Block", all.x = T)
  }else{
    spots$vshift <- 0
    spots$hshift <- 0
  }
  #apply shifts
  if(apply.shifts){
    spots <- rppa.vshift(spots)
    spots <- rppa.hshift(spots)
  }
  
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
  
  #filter bad spots (set signal to NA)
  if(filter.diameter) spots <- rppa.filter.diameter(spots)
  if(filter.flag) spots <- rppa.filter.flagged(spots)
  if(filter.neg.values) spots <- rppa.filter.neg.values(spots)
  
  blocksUrl <- paste(baseUrl, "getBlocksPerRow/", slideIndex, sep = "")
  if(!is.na(securityToken)) blocksUrl <- paste(blocksUrl, "?securityToken=", securityToken, sep="")
  
  titleUrl <- paste(baseUrl, "getTitle/", slideIndex, sep = "")
  if(!is.na(securityToken)) titleUrl <- paste(titleUrl, "?securityToken=", securityToken, sep="")
  
  antibodyUrl <- paste(baseUrl, "getAntibody/", slideIndex, sep = "")                     
  if(!is.na(securityToken)) antibodyUrl <- paste(antibodyUrl, "?securityToken=", securityToken, sep="")
  
  pmtUrl <- paste(baseUrl, "getPMT/", slideIndex, sep = "")                     
  if(!is.na(securityToken)) pmtUrl <- paste(pmtUrl, "?securityToken=", securityToken, sep="")
  
  spots <- rppa.set.blocksPerRow(spots, as.integer(scan(text=getURL(blocksUrl, curl=connection), what = "integer")))
  attr(spots, "PMT") <- as.integer(scan(text=getURL(pmtUrl, curl=connection), what = "integer"))
  spots <- rppa.set.title(spots, paste(scan(text=getURL(titleUrl, curl=connection), what = "character"), collapse=" "))
  spots <- rppa.set.antibody(spots, paste(scan(text=getURL(antibodyUrl, curl=connection), what = "character"), collapse=" "))
  attr(spots, "slideIndex") <- slideIndex
  message("...everything done. returning data.")
  return(spots)
}