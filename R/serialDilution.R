rppa.serialDilution.manipulate <- function(spots){
  require(manipulate)
  
  manipulate(
    View(rppa.serialDilution(spots, initial.dilution.estimate, sensible.min, sensible.max)),
    initial.dilution.estimate=slider(1, 5, step=0.1, initial=2),
    sensible.min = slider(1,1000, step=10, initial=1),
    sensible.max = slider(1000,100000, step=1000, initial=60000)
  )
}

rppa.serialDilution.batch <- function(slideList)
{
  require(plyr)
  
  for(slide in slideList)
  {
    if(is.null(attr(slide, "title"))){
      stop("One or more slides are without title! Please use rppa.set.title to assign a title before comparing multiple slides.")
      return()
    }
  }
  
  data.protein.conc <- ldply(slideList, function(x) { 
    result <- rppa.serialDilution(x)
    result$Slide <- attr(x, "title")
    return(result)
  })
  
  return(data.protein.conc)
}

rppa.serialDilution <- function(spots, initial.dilution.estimate=2, sensible.min=1, sensible.max=6e4, method="nls", compress.results=T, make.plot=T, useDepositions=F,...)
{ 
  #get title
  slideTitle <- attr(spots, "title")
  
  #subset to sample spots if possible
  sample.spots <- subset(spots, SpotClass=="Sample")
  if(!is.null(dim(sample.spots))) spots <- sample.spots
  
  #convert input table so that each dilution is in one column
  spots.c <- rppa.serialDilution.format(spots)
  
  #extract number of different dilutions that are not NA
  numOfDilutions <- length(unique(spots$DilutionFactor[!is.na(spots$DilutionFactor)])) 
  
  if(numOfDilutions <= 1) stop("Less than two dilutions found in data. Cannot compute a serial dilution curve.")
  
  #calculate matrix of dilutions
  spots.m <- rppa.serialDilution.dilutionMatrix(spots.c, numOfDilutions)
  
  depositionCorrection <- function(slide){
    deposFactors = sort(unique(slide$Deposition))
    result <- foreach(i=2:length(deposFactors), .combine=cbind) %do%
    {
      correction <- median(slide[slide$Deposition == deposFactors[i], "Signal"],na.rm=T) /
      median(slide[slide$Deposition == deposFactors[1], "Signal"],na.rm=T)
    }
    result <- as.vector(cbind(1, result))
    names(result) <- deposFactors
    
    return(result)
  }
  
  if(!is.null(spots.c$Deposition)) correctedDepositions <- spots.c$Deposition
  #correctedDepositions <- depositionCorrection(spots)[as.character(spots.c$Deposition)]
  
  #compute the actual protein estimates using the serial dilution method
  spots.e <- rppa.serialDilution.compute(spots.m, initial.dilution.estimate, sensible.min, sensible.max, method, slideTitle, make.plot, correctedDepositions)
  
  #fit 
  fit <- attr(spots.e, "fit")
  fittedData <- attr(spots.e, "fittedData") 
  pairedData <- attr(spots.e, "pairedData") 
  estDilFactor <- attr(spots.e, "estimatedDilutionFactor") 
  
  #combine estimates with signal information
  spots.result <- cbind(spots.c[,1:(ncol(spots.c)-numOfDilutions)], spots.e)
  
  #filter values under detection limit or saturated
  spots.result[!is.na(spots.result$xflag),"x.weighted.mean"] <- NA
  
  if(!compress.results){
    attr(spots.result, "fit") <- fit
    return(spots.result)
  } 
  
  spots.summarize <- rppa.serialDilution.summarize(spots.result, ...)

  spots.summarize$concentrations <- spots.summarize$x.weighted.mean
  spots.summarize$upper <- spots.summarize$x.err + spots.summarize$x.weighted.mean
  spots.summarize$lower <- spots.summarize$x.weighted.mean - spots.summarize$x.err
    
  spots.summarize <- spots.summarize[,!(colnames(spots.summarize) %in% c("x.weighted.mean", "x.err", "flags"))]
  attr(spots.summarize, "title") <- attr(spots, "title")
  attr(spots.summarize, "antibody") <- attr(spots, "antibody")
  attr(spots.summarize, "fit") <- fit
  attr(spots.summarize, "fittedData") <- fittedData
  attr(spots.summarize, "pairedData") <- pairedData
  attr(spots.summarize, "estimatedDilutionFactor") <- estDilFactor
  
  readout <- data.frame(concentrations=spots.summarize$readout,
                        upper=spots.summarize$readout.sem + spots.summarize$readout,
                        lower=spots.summarize$readout - spots.summarize$readout.sem)
  readout.centered <- readout / mean(readout$concentrations, na.rm=T)
  attr(spots.summarize, "readout") <- readout
  attr(spots.summarize, "readout.centered") <- readout.centered
  
  return(spots.summarize)
}

rppa.serialDilution.format <- function(spots, useDepositions=F) {  
  require(reshape)
  
  title <- attr(spots, "title")
  
  #make sure dilution factor and deposition is not NA
  if(length(unique(spots$DilutionFactor)) == 1 && is.na(unique(spots$DilutionFactor))) spots$DilutionFactor <- 1
  if(length(unique(spots$Deposition)) == 1 && is.na(unique(spots$Deposition))) spots$Deposition <- 1
  
  #filter NA values
  spots <- spots[!is.na(spots$DilutionFactor),]
  groupingCols <- setdiff(colnames(spots), c("Block", "id", "Row", "Column", "Signal", "surface", "BlockRow", "BlockColumn", "DilutionFactor", "FG", "BG", "Flag", "Diameter", "SGADesc", "SGBDesc", "SGCDesc", "hshift", "vshift"))
  
  #transform continuous into descrete
  spots$DilutionFactor <- as.factor(spots$DilutionFactor)
  #if(!useDepositions) groupingCols <- c(groupingCols, "Deposition")
  groupingCols <- paste(groupingCols, collapse = " + ")
  
  if(useDepositions) dilutions <- "DilutionFactor+Deposition"
  else dilutions <- "DilutionFactor"    
  
  spots.c <- cast(spots, as.formula(paste(groupingCols, dilutions, sep=" ~ ")), value="Signal", add.missing=FALSE, fun.aggregate="median", na.rm=T)
  
  if(useDepositions) attr(spots.c, "numOfDilutions") <- length(levels(spots$DilutionFactor)) * length(unique(spots$Deposition))  
  else attr(spots.c, "numOfDilutions") <- length(levels(spots$DilutionFactor))

  return(spots.c)
} 

rppa.serialDilution.pairColumns <- function(spots, startColumn=1, depositions=NA){
  
  pairedData <- data.frame(x=numeric(0), y=numeric(0), depos=numeric(0))
  
  for(i in startColumn:(ncol(spots)-1))
  {
    pairedData <- rbind(pairedData, cbind(x=spots[,i+1], y=spots[,i], depos=depositions))  
  }
  
  return(pairedData)
  
}

rppa.serialDilution.dilutionMatrix <- function(spots.c, numOfDilutions, highestDilutionFirst=T)
{  
  #extract dilution matrix for serial dilution curve algorithm
  spots.m <- as.matrix(spots.c[,(ncol(spots.c)-(numOfDilutions-1)):ncol(spots.c)] )
  
  #make sure order is correct
  if((mean(spots.m[,1], na.rm=T)< mean(spots.m[,2], na.rm=T) && highestDilutionFirst) || (mean(spots.m[,1], na.rm=T) > mean(spots.m[,2], na.rm=T) && !highestDilutionFirst))
  {
    spots.m <- spots.m[,ncol(spots.m):1]
  }

  return(spots.m)
}

rppa.serialDilution.filter <- function(data, sensible.min, sensible.max)
{
  #filter NA values
  data <- subset(data, !is.na(x) & !is.na(y))
  
  #filter bad log odds ratios
  ratio <- log(abs(data$y)/abs(data$x)) 
  ratio.median <- median(ratio, na.rm=T)
  ratio.mad <- mad(ratio,na.rm=T)
  filter <- abs(ratio - ratio.median)/2 > ratio.mad
  
  data <- data[!filter,] 
  
  #filter non-sensible values
  data <- subset(data, x > sensible.min & x < sensible.max & y > sensible.min & y < sensible.max)
  
  return(data)
}

rppa.serialDilution.compute <- function(spots.m, initial.dilution.estimate=2.5, sensible.min=5, sensible.max=6e4, method="nls", slideTitle="", make.plot=T, depositions=NA){

  #pair columns for serial dilution plot
  pairedData <- rppa.serialDilution.pairColumns(spots.m, depositions=depositions)
  
  #starting values for non-linear model fit
  a <- max(5, min(spots.m, na.rm=T))
  M <- max(spots.m, na.rm =T) 
  D0 <- initial.dilution.estimate
  D <- D0
  minimal.err <- 5
  
  #filter non-sensible values
  pairedData <- rppa.serialDilution.filter(pairedData, sensible.min, sensible.max)
  
  #nls model fit
  if(method=="nlsLM") fit <- nlsLM(y ~ a +1/((1/(x -a) -c)/(D)+c), data=pairedData, start=list(a=a,D=D,c=1/M))
  else fit <- nls(y ~ a +1/((1/(x -a) -c)/(D)+c), data=pairedData, start=list(a=a,D=D,c=1/M),alg="port", lower=list(minimal.err,1,0),weights=1/(minimal.err+abs(as.numeric(x))))    
  
  #calculate fitted data
  fittedData <- data.frame(x=pairedData$x, depos=pairedData$depos, y=predict(fit, data.frame(x=pairedData$x, depos=pairedData$depos)))
  #assemble parameters for serial dilution algorithm
  a <- summary(fit)$parameter[1]
  D <- summary(fit)$parameter[2]
  c <- summary(fit)$parameter[3]
  d.a <- summary(fit)$parameter[4]
  d.D <- summary(fit)$parameter[5]
  d.c <- summary(fit)$parameter[6]
  M <- a+1/summary(fit)$parameter[3]
  
  if(make.plot){
    #plot serial dilution curve
    library(ggplot2)
    library(gridExtra)
    print(ggplot(fittedData, aes(x=y, y=log2(x))) + labs(title=paste("Signal to concentration estimate plot: ", slideTitle)) + ylab("Signal") + xlab("Concentration estimate") + geom_point(aes(color=as.factor(depos))) + geom_smooth(aes(group=depos), method="loess"))
    sdc.plot <- ggplot(pairedData, aes(x=x, y=y)) + labs(title=paste("Serial Dilution Curve Fit: ", slideTitle, ", estimated dilution factor ", round(D, 2))) + xlab("Signal at next dilution step") + ylab("Signal") + geom_point(aes(color=as.factor(depos))) #+ coord_cartesian(ylim = c(0, sensible.max)) 
    sdc.plot <- sdc.plot + geom_abline(intercept=0, slope=1, color="red")

    sdc.plot <- sdc.plot + geom_line(data=fittedData)
    
    sdc.plot <- sdc.plot + scale_color_discrete(name="Depositions")    
    print(sdc.plot)
  }
  
  #estimate protein concentrations
  serialDilutionResult <- rppa.serialDilution.protein.con(D0=D0, D=D,c=c,a=a,d.a=d.a, d.D=d.D, d.c=d.c,data.dilutes=spots.m,depos=depositions)
  
  #add fit object
  attr(serialDilutionResult, "fit") <- fit
  attr(serialDilutionResult, "fittedData") <- fittedData
  attr(serialDilutionResult, "pairedData") <- pairedData
  attr(serialDilutionResult, "estimatedDilutionFactor") <- round(D,2)
  
  return(serialDilutionResult)
}

rppa.serialDilution.protein.con <- function (D0,D,c,a,d.D,d.c, d.a, data.dilutes,r=1.1,minimal.err=50, depos) {
  #D0 = dilution.factor # this is a preset value in diluteion experiments, typical D0=10, 3, or 2.
  #D fitted dilution factor. Ideally, D = D0 ^ gamma, where gamma is a parameter in Sips model
  # k 1:ncol(data.dilutes)    # index of dilute dilution steps in each dilution sereies
  # Np 1:nrow(data.dilutes)   # index of samples
  x.weighted.mean= rep(NA,nrow(data.dilutes))
  x.err = x.weighted.mean
  xflag = x.weighted.mean    # takes values of 0,NA,2, which means under detection, OK, saturated
  K = ncol(data.dilutes)   	# number of total dilution steps for a sample
  igamma = log(D0)/log(D)	#where gamma is 1/gamma, a parameter in Sips model
  M =min(6e4,1/c+a)			#when M is too large, take 1e9.
  
  x.saturation.level=   D0^(K-1)/((1/( M/r - a)- 1/(M-a)))^igamma 
  x.nodetection.level = D0^(1-1)/((1/( r*a - a)- 1/(M-a)))^igamma
  
  for (Np in 1:nrow(data.dilutes)){ # for each sample
    x=rep(NA,K); w=x; xL=x; xH=x;  #initialization
    y = data.dilutes[Np,]
    if(length(y[y<M/r])<2) {#condition to call saturation
      xflag[Np] = 2 
      x.weighted.mean[Np] = x.saturation.level # Use M/r value
      x.err[Np] = NA
    } else {
      if(length(y[y>r*a])<2) {#condition to call undetected
        
        xflag[Np] = 0
        x.weighted.mean[Np] = x.nodetection.level # Use r*a value
        x.err[Np] = NA
      } else {
        
        y[y>M/r] = NA # for removing signals near saturation 
        y[y<a*r] = NA # for removing signals near bg noise
        
        for (k in 1:K){# for each signal in a dilution series
          y[k] =max(min(M/1.01,y[k]), a+minimal.err) # limit y[k] to be within a+minimal.err and M/1.01
          x[k] =   D0^(k-1) * depos[Np] /(1/(y[k]-a)- c)^igamma #estimated protein concentration prior dilution
          #estimate the derivitives
          de.x.over.de.a = igamma * D0^(k-1) * depos[Np] *(1/(y[k]-a)- c)^(-igamma-1)/(y[k]-a)^2
          
          de.x.over.de.c = igamma * D0^(k-1) * depos[Np] * (1/(y[k]-a)- c)^(-igamma-1)
          de.x.over.de.D = x[k] *log(1/(y[k]-a)- c) * igamma/D/log(D)/(D0^(k-1) * depos[Np])
          w[k] = (de.x.over.de.a * d.a)^2 + ( de.x.over.de.c * d.c)^2 + (de.x.over.de.D * d.D)^2
        }
        w = w[!is.na(x)] # removing signals near saturation or bg noise
        x = x[!is.na(x)] # removing signals near saturation or bg noise
        
        if(length(x) > 0 ) {
          x.range = 3* max(1, median(x)*0.01,mad(x)) 
          x.f = (abs(x-median(x)) < x.range) # removing outliers
          #c(x,w)
          x=x[x.f]
          w=w[x.f] # removing outliers
          w= 1/w
          x.weighted.mean[Np] = sum (x*w) /sum(w)
          x.err[Np]=1/sqrt(sum(w))
        }
      }#end of else saturation
    }# end of else below detection 
    
  }#end of for each sample Np
  #return value:
  cbind(x.weighted.mean, x.err,xflag)
}#end of function

