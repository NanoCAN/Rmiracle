rppa.batch.sdc <- function(slideList, normalizeTo=NA, csv2=F, surface.normalization=T, plotHeatmaps=T, plotConcEstimates=T, saveDir=NA, positive.control="IgG 400",
                           swap=F, horizontal.line=T, fill.legend=T, error.bars=T, scales="free", sample.subset=NA, reference=NA,...)
{
  #requirements
  require(foreach)
  require(plyr)
  require(ggplot2)
  
  #processing serial dilution curve and surface normalization method
  sdc <- function(slide)
  {
    cat(paste("Processing", attr(slide, "title"), "...\n"))
    if(!is.na(saveDir)[1]){
      #heatmap
      if(plotHeatmaps){
        png(paste(attr(slide, "title"), "- Heatmap.png"), width=1024, height=768)
        rppa.plot.heatmap(slide, title=attr(slide, "title"))
        dev.off()
      }
      
      #write raw data
      if(csv2) write.csv2(slide, paste(attr(slide, "title"), "Raw Data.csv"))
      else write.csv(slide, paste(attr(slide, "title"), "Raw Data.csv"))
    }
  
    if(surface.normalization){
      cat("Applying surface normalization...\n")
      if(!is.na(saveDir)[1]) png(paste(attr(slide, "title"), "- Surface.png"), width=1024, height=768)  
      
      slide <- rppa.surface.normalization(slide, positive.control)
      
      if(!is.na(saveDir)[1]) dev.off()
      
      #heatmap after surface normalization
      if(plotHeatmaps && !is.na(saveDir)[1]){
        png(paste(attr(slide, "title"), "- Heatmap After Surface Normalization.png"), width=1024, height=768)
        rppa.plot.heatmap(slide, title=attr(slide, "title"))
        dev.off()
        
        png(paste(attr(slide, "title"), "- Heatmap Showing Surface Effects.png"), width=1024, height=768)
        rppa.plot.heatmap(slide, fill="surface", title=attr(slide, "title"))
        dev.off()
      }
    }
    
    #serial dilution curve plot
    if(!is.na(saveDir)[1]) png(paste(attr(slide, "title"), "- Serial Dilution Curve Fit.png"), width=1024, height=768)
    cat("Serial dilution curve fit...\n")
    
    #backup title
    slideTitle <- attr(slide, "title")
    slideAntibody <- attr(slide, "antibody")
    
    #call serial dilution curve
    slide <- rppa.serialDilution(slide, ...)
    
    #reset title
    attr(slide, "title") <- slideTitle
    attr(slide, "antibody") <- slideAntibody
  
    if(!is.na(saveDir)[1]){
      #close device for sdc plot
      dev.off()

      #write out estimates found through serial dilution curve fit
      if(csv2) write.csv2(slide, paste(attr(slide, "title"), "Serial Dilution Curve Estimates.csv"))
      else write.csv(slide, paste(attr(slide, "title"), "Serial Dilution Curve Estimates.csv"))
    }
    return(slide)
  }
  
  #check if data should be saved
  if(!is.na(saveDir)[1])
  {
    keepWd <- getwd()
    setwd(saveDir)
  }
  
  #method for protein concentration plot
  proteinConc <- function(pConc, normalized){
    if(!is.na(saveDir)[1]) png(paste(attr(pConc, "title"), "- Protein Concentration Estimates", normalized, ".png"), width=1024, height=768)
    rppa.proteinConc.plot(pConc, attr(pConc, "title"), swap, horizontal.line, fill.legend, error.bars, scales, sample.subset, reference)
    if(!is.na(saveDir)[1]) dev.off()
  }
  
  #serial dilution curve on all slides except normalization slide
  result.sdc <- foreach(slide=slideList) %do% sdc(slide)
  
  #plot protein concentration 
  if(plotConcEstimates)
  {
    cat("Plotting protein concentration estimates...\n")
    foreach(pConc=result.sdc) %do% proteinConc(pConc, "non-normalized")
    cat("...Done!\n")
  }
  
  #method for protein concentration normalization
  proteinConcNorm <- function(slide.sdc, normalizeTo.sdc){  
    slide.normalized <- rppa.proteinConc.normalize(slide.sdc, normalizeTo.sdc)
    slide.normalized$Slide <- attr(slide.sdc, "antibody")
    
    if(!is.na(saveDir)[1]){
      
      if(csv2) write.csv2(slide.normalized, paste(attr(slide.sdc, "title"), " - PConcEst normalized by", attr(normalizeTo.sdc, "antibody"), ".csv"))
      else write.csv(slide.normalized, paste(attr(slide.sdc, "title"), " - PConcEst normalized by", attr(normalizeTo.sdc, "antibody"), ".csv"))
    }
    return(slide.normalized)
  }
  
  #process normalization slide
  if(!is.na(normalizeTo)[1])
  {
    normalize <- function(normalizeTo, result.sdc){
      normalizeTo.sdc <- sdc(normalizeTo)
      cat(paste("Normalizing slides to", attr(normalizeTo, "title"), "...\n"))
      foreach(slide.sdc=result.sdc) %do% proteinConcNorm(slide.sdc, normalizeTo.sdc) 
    }
    
    if(!is.data.frame(normalizeTo) && is.list(normalizeTo))
    { 
      #TODO normalizeTo <- rppa.proteinConc.average(normalizeTo)
    }
    result.normalized <- normalize(normalizeTo, result.sdc)
  }
  
  if(exists("result.normalized")) result <- result.normalized
  else{
    result <- foreach(slide.sdc=result.sdc) %do%{
      slide.sdc$Slide <- attr(slide.sdc, "antibody")
      #center data 
      slide.sdc <- ddply(slide.sdc, .(A,B, Fill), function(slide.sdc){
        meanConcentrations <- mean(slide.sdc$concentrations, na.rm=T) 
        slide.sdc$concentrations <- slide.sdc$concentrations / meanConcentrations
        slide.sdc$upper <- slide.sdc$upper / meanConcentrations
        slide.sdc$lower <- slide.sdc$lower / meanConcentrations
        return(slide.sdc)
      })
      
      return(slide.sdc)
    }
  } 
  
  #plot protein concentration 
  if(plotConcEstimates && !is.na(normalizeTo[1]))
  {
    cat("Plotting protein concentration estimates (normalized)...\n")
    foreach(pConc=result) %do% proteinConc(pConc, paste("normalized to ", attr(normalizeTo, "antibody")))
    cat("...Done!\n")
  }
  if(!is.na(saveDir)[1]) setwd(keepWd)
  cat("Done!\n")
  
  cat("Plotting overview...\n")
  rppa.proteinConc.overview(ldply(result), title="Antibody Protein Estimate Comparison", subset.sample=sample.subset)
  cat("Everything done!\n")
  
  return(result)
}

rppa.batch.dunnett<- function(batch.result, referenceSample, p.cutoff=1, sample.subset=NA, duplicate.nas=T)
{
  #subset sample
  if(!is.na(sample.subset)[1]){
    batch.result <- lapply(batch.result, function(x, sample.subset, duplicate.nas){
      result <- subset(x, Sample %in% sample.subset)
      result <- rppa.duplicate.nas(result)
      result$Sample <- factor(result$Sample, sample.subset)
      return(result)
    }, sample.subset=sample.subset, duplicate.nas=duplicate.nas)
  }
  else{
    batch.result <- foreach(slide=batch.result) %do% rppa.duplicate.nas(slide)
  }
  
  #for each of the slides calculate the test statistics (pairwise comparison using t-test with multiple comparison adjustment)
  pvalues <- foreach(slide=batch.result, .combine=rbind) %dopar%
  {
    rppa.dunnett(slide, referenceSample)
  }
  pvalues$Samples <- factor(pvalues$Samples, paste(sample.subset, "-", referenceSample))
  
  rppa.batch.dunnett.plot(pvalues, p.cutoff)
  
  return(pvalues)
}

rppa.batch.dunnett.plot <- function(pvalues, p.cutoff=1)
{
  require(ggplot2)
  require(scales)
  pvalues.subset <- subset(pvalues, pvalues <= p.cutoff)
  pvalues.subset$symbol <- ""
  pvalues.subset[pvalues.subset$pvalues < 0.01, "symbol"] <- "*"
  pvalues.subset[pvalues.subset$pvalues < 0.001, "symbol"] <- "**"
  pvalues.subset[pvalues.subset$pvalues < 0.0001, "symbol"] <- "***"
  limits <- aes(ymax = estimates + stderror, ymin = estimates - stderror)
  q <- qplot(x=Samples, y=estimates, data=pvalues.subset, fill=pvalues, ylab="estimated difference", geom="bar", stat="identity", label=symbol)
  q <- q + geom_errorbar(limits, position="dodge", width=0.25)
  q <- q + theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1))
  q <- q + scale_fill_gradient2(trans="log", low="red", guide="legend", mid="orange", high="yellow", midpoint=1e-6, breaks=10^(-(seq(-3, 12, by=3))))
  q <- q + facet_grid(slide ~ A + B)
  q <- q + geom_text(aes(y = estimates + stderror), vjust=0.1)
  q <- q + scale_y_continuous(labels = percent)
  print(q)
}