rppa.tukeyHSD <- function(slide)
{
  foreach(A = levels(slide$A), .combine=rbind) %do%{
    foreach(B = levels(slide$B), .combine=rbind) %do%{
      slide.subset <- subset(slide, A==A & B==B)
      tukey.df <- as.data.frame(TukeyHSD(aov(concentrations ~ Sample, data=slide.subset), ordered=T)$Sample)
      tukey.df$Samples <- row.names(tukey.df)
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
  q <- q + scale_y_continuous(labels = percent)
  print(q)
}

rppa.dunnett <- function(slide, referenceSample="OTC#3")
{
  require(multcomp)
  foreach(currentA = levels(slide$A), .combine=rbind) %do%{
    foreach(currentB = levels(slide$B), .combine=rbind) %do%{
      slide.subset <- subset(slide, A==currentA & B==currentB)
      slide.subset$Sample <- relevel(slide.subset$Sample, ref=referenceSample)
      slide.aov <- aov(concentrations ~ Sample, data=slide.subset)
      dunnett.df <- with(summary(glht(slide.aov, linfct=mcp(Sample="Dunnett")))$test, 
                         {  data.frame(estimates=coefficients, stderror=sigma, pvalues=pvalues) } )
      dunnett.df$Samples <- row.names(dunnett.df)
      dunnett.df$A <- currentA
      dunnett.df$B <- currentB
      dunnett.df$slide <- slide$Slide[1]
      return(dunnett.df)
    }
  }
}

rppa.plot.dunnett <- function(pvalues, p.cutoff=1)
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
  q <- q + scale_fill_gradient2(trans="log", low="red", guide="legend", mid="orange", high="yellow", breaks=10^(-(seq(-3, 12, by=3))))
  q <- q + facet_grid(A ~ B)
  q <- q + geom_text(aes(y = estimates + stderror), vjust=0.1)
  q <- q + scale_y_continuous(labels = percent)
  print(q)
}