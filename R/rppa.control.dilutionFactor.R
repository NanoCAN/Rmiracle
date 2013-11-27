rppa.control.dilutionFactor <- function(spots, ncol=2, select.columns.A="CellLine", merge.plots=T, select.columns.B="Treatment", select.columns.fill="Deposition"){
  library(ggplot2)
  library(gridExtra)
  library(foreach)
  dilutions <- sort(unique(spots$DilutionFactor[!is.na(spots$DilutionFactor)]), decreasing=T)
  numOfDilutions <- length(dilutions) 
  test.c <- rppa.serialDilution.format(spots)
  test.m <- rppa.serialDilution.dilutionMatrix(test.c, numOfDilutions=numOfDilutions,highestDilutionFirst=T)
  test.e <- test.c[,1:(ncol(test.c)-numOfDilutions)]
  listOfRatios <- list()
  listOfPlots <- foreach(i=1:(ncol(test.m)-1)) %do%
{
  title <- paste("Comparing dilution ", dilutions[i], " and ", dilutions[i+1])
  ratio <- (test.m[,i] / test.m[,i+1])
  cat(title)
  cat("\n")
  print(summary(ratio))
  listOfRatios[[title]] <-  ratio
  temp.result <- test.e
  temp.result$ratio <- ratio
  q <- qplot(x=Deposition, y=ratio, data=temp.result, fill=as.factor(Deposition), main=title, geom="boxplot") 
  q <- q +  guides(fill = guide_legend(title = select.columns.fill))
  q <- q + geom_smooth(aes(group=1,fill = factor(Deposition)), method="loess")
  q <- q + facet_grid(CellLine ~ Treatment, scales="free")
  q <- q + stat_summary(aes(label=round(..y..,1)), fun.y=median, geom="text", size=6, vjust=-0.5)
  return(q)
}
  if(merge.plots) do.call("grid.arrange", c(listOfPlots, ncol=ncol, main=paste("Boxplots of signal increase between different dilutions in", attr(spots, "title"))))
  else lapply(listOfPlots, function(p) print(p))
  return(listOfRatios)
}