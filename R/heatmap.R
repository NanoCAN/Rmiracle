rppa.plot.heatmap <- function(spots, log=NA, fill="Signal", plotNA=T, palette="Set1", 
                              discreteColorA=NA, discreteColorB=NA, discreteColorC=NA, title=NA){
  
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
  
  if(nrow(spots[!is.na(spots[,fill]),]) == 0){
    stop("There is no information on that property!")
  } 
  
  if(is.na(title)){
   if(is.null(attr(spots, "title"))) title <- "" 
   else title <- attr(spots, "title") 
  }
  #fix row and column nums
  spots$Row <- as.integer(spots$Row)
  spots$Column <- as.integer(spots$Column)
  
  #transform continuous into descrete
  spots$Deposition <- as.factor(spots$Deposition)
  spots$Dilution <- as.factor(spots$Dilution)
  
  #reverse order of levels so that first factor will be purest concentration 
  spots$Dilution <- factor(spots$Dilution, levels=rev(levels(spots$Dilution)))
  
  if(!is.na(log))
  {
    if(log=="log2") spots$Signal <- log2(spots[[fill]])  
    else if(log=="log10") spots$Signal <- log10(spots[[fill]])
  }
  else{
    spots$Signal <- spots[[fill]]
  }
  
  p <- ggplot(spots, aes(x=Column, y=Row)) + ggtitle(title)
  
  if(plotNA) 
  {
    p <- p + geom_raster(data=subset(spots, !is.na(Signal)), aes_string(fill = fill));
    if(nrow(subset(spots,is.na(Signal)))!=0)
    {
      p <- p + geom_raster(data=subset(spots, is.na(Signal)), fill="black");
    }
  }
  else {
    p <- p + geom_raster(data=spots, aes_string(fill = fill));
  }
  
  #how many blocks per row?
  blocksPerRow <- 12
  if(!is.null(attr(spots, "blocksPerRow"))) {
    blocksPerRow <- attr(spots, "blocksPerRow") 
  }
  
  p <- p + coord_cartesian(ylim=c(max(spots$Row)+0.5,0.5));
  p <- p + facet_wrap(~Block, ncol=blocksPerRow);
  p <- p + scale_x_continuous(expand=c(0,0), breaks=seq(1, max(spots$Column), 3)) + scale_y_continuous(expand=c(0,0), trans="reverse");
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.margin=unit(0.1, "lines"), panel.margin=unit(0, "lines"), 
                plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"),   
                plot.title=element_text(size=18), strip.background=element_rect(fill="grey90", colour="grey50"))
  
  #colorbrewer color palette
  if(!is.na(palette) && is.factor(spots[[fill]])){
    getPalette <- colorRampPalette(brewer.pal(9, palette))
    discreteFactor <- spots[[fill]]
    p <- p + scale_fill_manual(values = getPalette(length(levels(discreteFactor))))
  }  
  
  else if(!is.na(discreteColorA) && !is.na(discreteColorB) && !is.factor(spots[[fill]]))
  {
    if(!is.na(discreteColorC)){
      p <- p + scale_fill_gradient2(low = discreteColorA, high = discreteColorB, mid=discreteColorC);
    }
    else
    {
      p <- p + scale_fill_gradient(low = discreteColorA, high = discreteColorB);
    }
  }
  
  print(p);
}


rppa.heatmap <- function(spots){  
  
  require(manipulate)
  
  manipulate(
  rppa.plot.heatmap(spots, log, fill, plotNA, palette, discreteColorA, discreteColorB, NA, NA)
, log = picker("none", "log2", "log10")
, fill = picker("Signal", "FG", "BG","Deposition", "CellLine", "LysisBuffer", "DilutionFactor", "Inducer", "SpotType", "SpotClass", "SampleName", "SampleType", "TargetGene")
, palette = picker("Set1", "Set2", "Set3", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2")
, plotNA = checkbox(TRUE, "Plot NA values")
, discreteColorA = picker("darkblue", "red", "blue", "steelblue", "magenta", "yellow", "white", "green")
, discreteColorB = picker("red", "darkblue", "blue", "steelblue", "magenta", "yellow", "white", "green")
  )
}


