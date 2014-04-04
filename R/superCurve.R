rppa.superCurve.parse.data <- function(Sample, spots)
{
  blocksPerRow <- attr(spots, "blocksPerRow")
  Sample <- apply(Sample, 1, paste, collapse=" # ")
  Sample[spots$SpotType!="Sample"] <- "Control"
  Mean.Net <- spots$Signal
  Mean.Total <- spots$FG
  Vol.Bkg <- spots$BG
  Main.Row <- (spots$Block %/% blocksPerRow) + 1
  Main.Col <- (spots$Block %% blocksPerRow) 
  
  #columns zero are actually columns blocksPerRow
  Main.Row[Main.Col==0] <- Main.Row[Main.Col==0] - 1
  Main.Col[Main.Col==0] <- blocksPerRow  
  Sub.Row <- spots$Row
  Sub.Col <- spots$Column
  
  parsed.data <- data.frame(Main.Row, Main.Col, Sub.Row, Sub.Col, Sample, Mean.Net, Mean.Total, Vol.Bkg)
  
  return(parsed.data)
}

rppa.superCurve.create.rppa <- function(parsed.data, spots)
{
  new.rppa = new("RPPA")
  new.rppa@data <- parsed.data
  new.rppa@file <- attr(spots, "title")
  new.rppa@antibody <- attr(spots, "antibody")
  
  return(new.rppa)
}

rppa.superCurve.create.df <- function(new.fit, groupingCols, log2=F)
{
  new.fit <- data.frame(x.weighted.mean=2^new.fit@concentrations, x.err=(2^new.fit@upper - 2^new.fit@concentrations), xflag=NA)
  
  new.cols <- strsplit2(row.names(new.fit), " # ")
  new.cols <- as.data.frame(new.cols)
  colnames(new.cols) <- groupingCols
  
  #fix NAs
  new.cols[new.cols=="NA"] <- NA
  new.cols <- apply(new.cols, 2, factor)
  
  new.df <- cbind(new.cols, new.fit)
  row.names(new.df) <- NULL
  return(new.df)
}

rppa.superCurve <- function(spots, return.fit.only=F, model="logistic", 
                            method="nls", ci=F, use.depositions=F, make.plot=T, trim=2, verbose=F, select.columns.A, select.columns.B, select.columns.fill){   
  library(limma)
  library(SuperCurve)
  
  #check for necessary attributes title, antibody, 
  if(is.null(attr(spots, "title"))) stop("Please set attribute 'title' first!")
  if(is.null(attr(spots, "antibody"))) stop("Please set attribute 'antibody' first!")
  if(is.null(attr(spots, "blocksPerRow")))stop("Please set attribute 'blocksPerRow' first!")
  
  #correct dilution factors
  spots$DilutionFactor <- as.double(spots$DilutionFactor)
  spots$Signal <- spots$FG #performs its own background estimation and can't deal so well with too many NAs
  
  #create data object for SuperCurve package  
  
  #create sample name from all selected columns  
  groupingCols <- setdiff(colnames(spots), c("Block", "id", "Row", "Column", "Signal", "surface", "BlockRow", "BlockColumn", "DilutionFactor", "FG", "BG", "Flag", "Diameter", "SGADesc", "SGBDesc", "SGCDesc", "hshift", "vshift"))
  if(!use.depositions) groupingCols <- c(groupingCols, "Deposition")
  Sample <- spots[,groupingCols]
  parsedData <- rppa.superCurve.parse.data(Sample, spots)
  
  #put the information in a RPPA data object
  new.rppa <- rppa.superCurve.create.rppa(parsedData, spots)
  
  #we need the dilution factors as log2 to a reference point, which we will choose to be undiluted 1.0
  steps <- round(log2(spots$DilutionFactor),2)
  if(use.depositions) steps <- steps + log2(spots$Deposition)
  steps <- steps 
  steps[is.na(steps)] <- 0
    
  new.design <- RPPADesign(new.rppa, steps=steps, series=new.rppa@data$Sample, controls=list("Control"), center=F, ordering="increasing")
  
  if(make.plot) image(new.design)
  
  new.fit <- RPPAFit(new.rppa, new.design, "Mean.Total", ci=ci, method=method, model=model, trim=trim,verbose=verbose)
  if(return.fit.only) return(new.fit)
  
  if(make.plot){
    plot(new.fit)
    image(new.fit)
  }

  new.df <- rppa.superCurve.create.df(new.fit, groupingCols)
  attr(new.df, "title") <- attr(spots, "title")
  attr(new.df, "antibody") <- attr(spots, "antibody")
  new.df <- new.df[new.df$SpotType%in%c("Sample"),]
  spots.summarize <- rppa.serialDilution.summarize(new.df, select.columns.A=select.columns.A, select.columns.B=select.columns.B, select.columns.fill=select.columns.fill)
  spots.summarize$concentrations <- spots.summarize$x.weighted.mean
  spots.summarize$upper <- spots.summarize$x.err + spots.summarize$x.weighted.mean
  spots.summarize$lower <- spots.summarize$x.weighted.mean - spots.summarize$x.err
  spots.summarize <- spots.summarize[!is.na(spots.summarize$Sample),]
  attr(spots.summarize, "fit") <- new.fit
  return(spots.summarize)
}