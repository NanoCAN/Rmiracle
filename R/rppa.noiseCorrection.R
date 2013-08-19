rppa.noiseCorrection <- function(dilutionMatrix, highestDilutionFirst=F, subtractNoiseFromAll=T){

  noisySpots <- data.frame()
  
  columnIterator <- function(colIndex, nextColIndex){
    if(!is.na(dilutionMatrix[row, colIndex]) && !is.na(dilutionMatrix[row, nextColIndex]))
    {
      if(dilutionMatrix[row,col] >= dilutionMatrix[row,nextColIndex]){
        return(TRUE)
      } 
    }
    return(FALSE)
  }
  
  for(row in 1:nrow(dilutionMatrix))
  {
    noisySpot <- FALSE
    if(!highestDilutionFirst){
      for(col in 1:(ncol(dilutionMatrix)-1))
      {
        if(columnIterator(col, (col+1)))
        {
          noisySpots<-rbind(noisySpots, c(row, col, dilutionMatrix[row, col]))
          dilutionMatrix[row,col] <- NA
        }
      }
    }
    else{
      for(col in ncol(dilutionMatrix):2)
      {
        if(columnIterator(col, (col-1)))
        {
          noisySpots<-rbind(noisySpots, c(row, col, dilutionMatrix[row, col]))
          dilutionMatrix[row,col] <- NA
        }
      }
    }
  }

  if(nrow(noisySpots)==0){
    cat("Could not detect noise.\n")
    return(dilutionMatrix)
  }
  noiseLevel <- summary(noisySpots[,3])[5]
  cat(paste(nrow(noisySpots), "spots detected where the higher diluted sample had a stronger signal.\n"))
  
  cat(paste("Noise detected:", noiseLevel, "\n"))
  if(subtractNoiseFromAll)
  {
    cat("Correcting signal values for noise level...\n")
    dilutionMatrix <- dilutionMatrix - noiseLevel
    dilutionMatrix[dilutionMatrix <= 0] <- NA
  }
  
  return(dilutionMatrix)
}