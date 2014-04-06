rppa.reformatColTypes <- function(spots)
{  
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
  spots$Replicate <- as.integer(spots$Replicate)
  spots$SpotType <- as.factor(spots$SpotType)
  spots$SpotClass <- as.factor(spots$SpotClass)
  
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

rppa.loadFromFile <- function(fileName=NULL, filter.diameter=T, filter.neg.values=T, filter.flag=T){
 
  #load files
  if(is.null(fileName)) fileName <- file.choose()
  spots <- read.delim(fileName, skip=8,na.string=c("NA", "null"))
  
  #apply shifts
  spots <- rppa.vshift(spots)
  spots <- rppa.hshift(spots)
  
  #skip two lines
  con <- file(fileName, open="r")
  readLines(con,2)
  
  #apply depositions
  depositionPattern <- readLines(con,1)
  depositionPattern <- gsub("\\[", "", depositionPattern)
  depositionPattern <- gsub("\\]", "", depositionPattern)
  depositionPattern <- as.integer(strsplit(depositionPattern, 
                                           ",")[[1]])
  spots$Deposition <- spots$Column%%length(depositionPattern)
  spots$Deposition[spots$Deposition == 0] <- length(depositionPattern)
  spots$Deposition <- depositionPattern[spots$Deposition]

  #format columns
  spots <- rppa.reformatColTypes(spots)
  
  #filter bad spots (set signal to NA)
  if(filter.diameter) spots <- rppa.filter.diameter(spots)
  if(filter.flag) spots <- rppa.filter.flagged(spots)
  if(filter.neg.values) spots <- rppa.filter.neg.values(spots)
  
  attr(spots, "slideIndex") <- readLines(con,1)
  spots <- rppa.set.title(spots, readLines(con,1))
  spots <- rppa.set.antibody(spots, readLines(con,1))
  attr(spots, "PMT") <- as.integer(readLines(con,1))
  spots <- rppa.set.blocksPerRow(spots, as.integer(readLines(con,1)))
  
  message("...everything done. returning data.")
  close(con)
  
  return(spots)
}

