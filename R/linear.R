rppa.linear <- function(slide, select.columns.sample="SampleName", select.columns.A, select.columns.B, select.columns.fill){
  library(stringr)
  slide <- subset(slide, SpotClass=="Sample")
  slide$DilutionFactor <- as.factor(slide$DilutionFactor)
  slide$Deposition <- as.factor(slide$Deposition)

  if(length(nchar(select.columns.sample)) == 1) selectedSamples <- as.character(slide[,select.columns.sample])
  else selectedSamples <- apply(slide[,select.columns.sample],1, paste, collapse=" | ")
  
  if(length(nchar(select.columns.A)) == 1) selA <- as.character(slide[,select.columns.A])
  else selA <- apply(slide[,select.columns.A],1, paste, collapse=" | ")
  
  if(length(nchar(select.columns.B)) == 1) selB <- as.character(slide[,select.columns.B])
  else selB <- apply(slide[,select.columns.B],1, paste, collapse=" | ")              
  
  if(length(nchar(select.columns.fill)) == 1) selFill <- as.character(slide[,select.columns.fill])
  selFill <- apply(slide[,select.columns.fill],1, paste, collapse=" | ")
  
  Sample <- apply(cbind(as.character(slide$SampleName), selA, selB, selFill), 1, paste, collapse="¤")
  df <- cbind(slide[,c("Signal", "DilutionFactor")], Sample)
  slide.rlm <- rlm(log(Signal) ~ Sample + DilutionFactor, data=df, maxit=100)
  samples <- grepl("Sample", names(slide.rlm$coefficients))
  result <- coef(summary(slide.rlm))[samples,]
  if(min(result[,1], na.rm=T) < 0)
    concentrations <- result[,1] - min(result[,1], na.rm=T)
  
  result <- data.frame(concentrations=concentrations, upper=concentrations + concentrations * result[,2], lower=concentrations - concentrations * result[,2])
  recreateSamples <- as.data.frame(t(matrix(unlist(str_split(row.names(result),pattern="¤")),4)))
  colnames(recreateSamples) <- c("Sample", "A", "B", "Fill")
  result <- cbind(recreateSamples, result)
  result[,"Sample"] <- sub("Sample", "", result[,"Sample"])
  return(result)
}

