rppa.serialDilution.summarize <-
function(data.protein.conc, method="mean", filter.bad.results=T,
                                          select.columns.sample=c("SampleName"), 
                                          select.columns.A=c("Treatment"),
                                          select.columns.B=c("CellLine"), select.columns.fill=c("NumberOfCellsSeeded"), useDeposition=T)
{
  data.protein.conc$x.err.percent <- data.protein.conc$x.err / data.protein.conc$x.weighted.mean
  
  
  if(length(select.columns.sample) > 1) Sample <- data.frame(Sample=apply(data.protein.conc[,select.columns.sample], 1 , paste, collapse=" | "))
  else Sample <- data.frame(Sample=data.protein.conc[,select.columns.sample])
  
  newColumns <- Sample
  
  if(!is.na(select.columns.A))
  {
    if(length(select.columns.A) > 1) A <- apply(data.protein.conc[,select.columns.A], 1 , paste, collapse=" | ")
    else A <- data.protein.conc[,select.columns.A]
    
    newColumns = cbind(newColumns, A)
  }
  
  if(!is.na(select.columns.B))
  {
    if(length(select.columns.B) > 1) B <- apply(data.protein.conc[,select.columns.B], 1 , paste, collapse=" | ")
    else B <- data.protein.conc[,select.columns.B]
    
    newColumns = cbind(newColumns, B)
  }
  
  if(!is.na(select.columns.fill))
  {
    if(length(select.columns.fill) > 1) Fill <- apply(data.protein.conc[,select.columns.fill], 1, paste, collapse=" | ")
    else Fill <- data.protein.conc[,select.columns.fill]
    
    newColumns = cbind(newColumns, Fill)
    fillAttribute <- paste(select.columns.fill, collapse=" | ")
  }
  else{
    cat("Fill attribute cannot be empty in the current version, using Replicate per default.")
    Fill <- data.protein.conc[,"Replicate"]
    newColumns <- cbind(newColumns, Fill)
    fillAttribute <- "Replicate"
  }

  if(filter.bad.results){
    data.protein.conc[!is.na(data.protein.conc$xflag), "x.weighted.mean"] <- NA
  } 
  selection <- c("x.weighted.mean", "x.err.percent", "PlateReadout")
  if(!("PlateReadout" %in% colnames(data.protein.conc))) data.protein.conc$PlateReadout <- NA
  if(useDeposition) selection <- c(selection, "Deposition")
  result <- cbind(newColumns, data.protein.conc[,selection])
  
  selection <- c("Sample", "Fill")
  if(!is.na(select.columns.A)) selection <- c(selection, "A")
  if(!is.na(select.columns.B)) selection <- c(selection, "B")
  if(useDeposition) selection <- c(selection, "Deposition")
  
  result <- ddply(result, selection, summarise, readout=mean(PlateReadout, na.rm=T), readout.sem=sqrt(var(PlateReadout,na.rm=TRUE)/length(na.omit(PlateReadout))), x.weighted.mean=mean(x.weighted.mean, na.rm=T), x.err=sum(x.err.percent, na.rm=T), sem=sqrt(var(x.weighted.mean,na.rm=TRUE)/length(na.omit(x.weighted.mean))))
  result$x.err <- result$x.weighted.mean * result$x.err
  
  #make sure A, B and sample are factors
  result$Sample <- as.factor(result$Sample)
  if(!is.na(select.columns.A)) result$A <- as.factor(result$A)
  if(!is.na(select.columns.B)) result$B <- as.factor(result$B)
  result$Fill <- as.factor(result$Fill)

  return(result)
}
