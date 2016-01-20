rppa.proteinConc.plot <- function(data.protein.conc, title="", swap=F, horizontal.line=T, fill.legend=T, error.bars=T, scales="free", sample.subset=NULL, reference=NULL, slideAsFill=F, ...){
  
  library(ggplot2)
  library(gridExtra) 
  library(dplyr)
  
  #subset samples and reorder 
  if(!is.null(sample.subset))
  {
    data.protein.conc <- subset(data.protein.conc, Sample %in% sample.subset)
    data.protein.conc$Sample <- factor(data.protein.conc$Sample, sample.subset)
  }
  
  #normalize data
  if(!is.null(reference)){
    data.protein.conc <- rppa.normalize.to.ref.sample(data.protein.conc, reference,...)    
  }
  
  #plot protein concentrations  
  limits <- aes(ymax = upper, ymin= lower)
  dodge <- position_dodge(width=0.9)
  if(slideAsFill) fill <- "Slide"
  else fill <- "Fill"

  p <- ggplot(data.protein.conc, aes_string(x="Sample", y="concentrations", fill=fill)) + ggtitle(title) +
                geom_bar(stat="identity", position="dodge") +
                ylab("Estimated Protein Concentration (Relative Scale)") +
                xlab("Sample")
              
  if(!is.null(data.protein.conc$Deposition)){  
    data.protein.conc <- rppa.mean.depos(data.protein.conc)
  }
  else if(!is.null(data.protein.conc$Fill))
  {
    data.protein.conc$Fill <- as.factor(data.protein.conc$Fill)
    if(fill.legend){
      p <- p + guides(fill=guide_legend(title=NULL))
    }
    else{
      p <- p + guides(fill=FALSE)
    }
  }
  else { 
    fill <- NULL
  }
  
  if(!is.null(data.protein.conc$B) && !is.null(data.protein.conc$A)){
    if(swap) p <- p + facet_grid(B~A, scales=scales)
    else p <- p + facet_grid(A~B, scales=scales)
  }
  
  else if(!is.null(data.protein.conc$A))
    p <- p + facet_wrap(~A, scales=scales)
  
  else if(!is.null(data.protein.conc$B))
    p <- p + facet_wrap(~B, scales=scales)
  
  if(error.bars) p <- p + geom_errorbar(limits, width=0.25, position=dodge)
  
  if(horizontal.line)  p <- p + geom_hline(aes(yintercept=1))
  
  p <- p + theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1))
  p <- p + theme(plot.margin = unit(c(1,2,1,1), "cm"))
  
  print(p)
  return(data.protein.conc)
}

rppa.mean.depos <- function(data.protein.conc){
  library(dplyr)
  data.protein.conc$Deposition <- as.factor(data.protein.conc$Deposition)
  #convert error bars to percentage
  data.protein.conc <- mutate(data.protein.conc, upper=(upper-concentrations)/concentrations, lower=(concentrations-lower)/concentrations)
  data.protein.conc <- data.protein.conc %>% regroup(lapply(intersect(colnames(data.protein.conc), 
                                                                      c("Sample", "Fill", "A", "B")), as.symbol)) %>% 
    summarise(concentrations=mean(concentrations, na.rm=T), 
              upper=sum(upper, na.rm=T), lower=sum(lower, na.rm=T))
  #revert to absolute errors
  data.protein.conc <- mutate(data.protein.conc, upper=(1+upper)*concentrations, lower=(1-lower)*concentrations)
  return(data.protein.conc)
}