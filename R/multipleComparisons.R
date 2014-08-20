rppa.tukeyHSD <- function(slide, control="OTC#3")
{
  levelsA = levels(as.factor(as.character(slide$A)))
  if(is.null(levels(slide$A))) levelsA <- "NA"
  
  levelsB = levels(as.factor(as.character(slide$B)))
  if(is.null(levels(slide$B))) levelsB <- "NA"
  
  foreach(A = levels(slide$A), .combine=rbind) %do%{
    foreach(B = levels(slide$B), .combine=rbind) %do%{
      slide.subset <- subset(slide, A==A & B==B)
      tukey.df <- as.data.frame(TukeyHSD(aov(Signal ~ SampleName, data=slide.subset), ordered=T)$SampleName)
      
      tukey.df$Samples <- row.names(tukey.df)
      tukey.df <- tukey.df[grepl(control, tukey.df$Samples),]
      
      avgControl <- mean(subset(slide.subset, SampleName %in% control)$Signal, na.rm=T)
      tukey.df$diff <- tukey.df$diff + avgControl
      tukey.df$upr <- tukey.df$upr + avgControl
      tukey.df$lwr <- tukey.df$lwr + avgControl
      
      tukey.df$A <- A
      tukey.df$B <- B
      tukey.df$slide <- slide$Slide[1]
      return(tukey.df)
    }
  }
}

rppa.plot.tukey <- function(pvalues, p.cutoff=1)
{
  require(ggplot2)
  require(scales)
  pvalues.subset <- subset(pvalues, `p adj` <= p.cutoff)
  pvalues.subset$symbol <- ""
  pvalues.subset[pvalues.subset[["p adj"]] < 0.01, "symbol"] <- "*"
  pvalues.subset[pvalues.subset[["p adj"]] < 0.001, "symbol"] <- "**"
  pvalues.subset[pvalues.subset[["p adj"]] < 0.0001, "symbol"] <- "***"
  limits <- aes(ymax = upr, ymin = lwr)
  q <- qplot(x=Samples, y=diff, data=pvalues.subset, fill=`p adj`, ylab="estimated difference", geom="bar", stat="identity", label=symbol)
  q <- q + geom_errorbar(limits, position="dodge", width=0.25)
  q <- q + theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1))
  q <- q + scale_fill_gradient2(trans="log", low="red", guide="legend", mid="orange", high="yellow", breaks=10^(-(seq(-3, 12, by=3))))
  q <- q + facet_grid(A ~ B)
  #q <- q + geom_text(aes(y = diff + stderror), vjust=0.1)
  #q <- q + scale_y_continuous(labels = percent)
  print(q)
}

rppa.dunnett <- function(slide, referenceSample="OTC#3", sample.subset=NULL)
{
  if(!is.null(sample.subset))
  slide <- subset(slide, Sample %in% sample.subset)
  
  levelsA = levels(as.factor(as.character(slide$A)))
  if(is.null(levels(slide$A))) levelsA <- "NA"

  levelsB = levels(as.factor(as.character(slide$B)))
  if(is.null(levels(slide$B))) levelsB <- "NA"
    
  library(multcomp)
  library(foreach)
  
  foreach(currentA = levelsA, .combine=rbind) %do%{
    foreach(currentB = levelsB, .combine=rbind) %do%{
      if(currentA == "NA" && currentB == "NA") slide.subset <- slide
      else if(currentA == "NA") slide.subset <- subset(slide, B==currentB)
      else if(currentB == "NA") slide.subset <- subset(slide, A==currentA)
      else slide.subset <- subset(slide, A==currentA & B==currentB)
      #check if ref exists
      if(nrow(subset(slide.subset, Sample == referenceSample)) == 0) stop(paste("Reference sample could not be found in group", currentA, "/", currentB, ". Change grouping parameters or select a subset of samples."))
      
      #center data first
      #slide.subset$concentrations <- slide.subset$concentrations / mean(slide.subset$concentrations, na.rm=T)
      
      slide.subset$Sample <- relevel(as.factor(as.character(slide.subset$Sample)), ref=referenceSample)
      slide.aov <- aov(concentrations ~ Sample, data=slide.subset)
   
      dunnett.df <- with(summary(glht(slide.aov, linfct=mcp(Sample="Dunnett")))$test, 
                         {  data.frame(estimates=coefficients, stderror=sigma, pvalues=pvalues) } )
      
      dunnett.df$Samples <- row.names(dunnett.df)
      if(!(currentA == "NA")) dunnett.df$A <- currentA
      if(!(currentB == "NA")) dunnett.df$B <- currentB
      if(!is.null(slide$Slide)) dunnett.df$slide <- slide$Slide[1]
      return(dunnett.df)
    }
  }
}

rppa.plot.dunnett <- function(pvalues, p.cutoff=1, set.neg.control.to.one=F, title="Dunnett's test")
{
  require(ggplot2)
  require(scales)
  
  pvalues.subset <- subset(pvalues, pvalues <= p.cutoff)  
  
  if(set.neg.control.to.one) pvalues.subset$estimates <- (pvalues.subset$estimates + 1)
  
  pvalues.subset$symbol <- ""
  pvalues.subset[pvalues.subset$pvalues < 0.05, "symbol"] <- "*"
  pvalues.subset[pvalues.subset$pvalues < 0.01, "symbol"] <- "**"
  pvalues.subset[pvalues.subset$pvalues < 0.001, "symbol"] <- "***"
  
  breaks <- c(c(1, 0.1, 0.05, 0.01, 0.001), 10^(-(seq(6, 21, by=3))))
  limits <- aes(ymax = estimates + stderror, ymin = estimates - stderror)
  q <- qplot(x=Samples, y=estimates, data=pvalues.subset, main=title, fill=pvalues, ylab="estimated difference", geom="bar", stat="identity", label=symbol)
  q <- q + geom_errorbar(limits, position="dodge", width=0.25)
  q <- q + theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1))
  q <- q + scale_fill_gradient2(low="red", trans="log", mid="orange", high="yellow", breaks=breaks, na.value="darkred", guide=guide_legend(keyheight=1), labels=paste(breaks, c("", "", ".", "*", "**", rep("***",6))))
  if(!is.null(pvalues$A) && !is.null(pvalues$B)) q <- q + facet_grid(A ~ B)
  else if(is.null(pvalues$A) && !is.null(pvalues$B)) q <- q + facet_grid(~B)
  else if(is.null(pvalues$B) && !is.null(pvalues$A)) q <- q + facet_grid(~A)
  q <- q + geom_text(aes(y = estimates + stderror), vjust=0.1)
  q <- q + scale_y_continuous(labels = percent)
  print(q)
}