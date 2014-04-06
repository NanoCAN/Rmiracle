library(shiny)
library(shinyIncubator)
library(stringr)
library(Rmiracle)
library(reshape2)

source("source//multiplot.R")

shinyServer(function(input, output, session) {
  
  ### DATA INPUT ###
  
  # Parse the GET query string
  output$queryText <- renderText({
    query <- parseQueryString(session$clientData$url_search)
    
    # Return a string with key-value pairs
    paste(names(query), query, sep = "=", collapse=", ")
  })
  
  getSlide <- function(baseUrl, securityToken){
    rppa.load(baseUrl=baseUrl, securityToken=securityToken)
  }  
  
  loadAllSlides <- reactive({
    
    query <- parseQueryString(session$clientData$url_search)
    if(length(query$baseUrl > 0)) baseUrl <- query$baseUrl
    else baseUrl <- "http://localhost:8080/MIRACLE/spotExport/"
    if(query$slideSecurityTokens != "") slideTokens <- str_split(query$slideSecurityTokens, "\\|")[[1]]
    else return(NULL)
    
    withProgress(session, min=1, max=(length(slideTokens)+1), expr={
      setProgress(message = 'Fetching readout data from MIRACLE...',
                  detail = 'Please be patient!',
                  value=0)
      readouts <- loadAllReadouts()
      slides <- list()
      count <- 1
      for(token in slideTokens)
      {
        setProgress(message = 'Fetching slide data from MIRACLE...',
                  detail = 'Please be patient!',
                  value=count)
        currentSlide <- getSlide(baseUrl, token)
        
        newSlide <- currentSlide
        if(!is.null(readouts)){
          newSlide <- merge(currentSlide, readouts[,c("PlateLayout", "PlateRow", "PlateCol", "PlateReadout")], all.x=T, by=c("PlateLayout", "PlateRow", "PlateCol"))
          #attributes(newSlide) <- attributes(currentSlide)
          attr(newSlide, "antibody") <- attr(currentSlide, "antibody")
          attr(newSlide, "slideIndex") <- attr(currentSlide, "slideIndex")
          attr(newSlide, "title") <- attr(currentSlide, "title")
          attr(newSlide, "blocksPerRow") <- attr(currentSlide, "blocksPerRow")
          #colnames(newSlide)[length(colnames(newSlide))] <- "PlateReadout"
        } 
        if(nrow(currentSlide) != nrow(newSlide)) stop("Could not match readout data. Too many matches")
        else{
          currentSlide <- newSlide
        } 
        slides[[attr(currentSlide, "slideIndex")]] <- currentSlide
        count <- count + 1
      }
      return(slides)
    })
  })
  
  loadAllReadouts <- reactive({  
    query <- parseQueryString(session$clientData$url_search)
    if(length(query$baseUrl > 0)) baseUrl <- query$baseUrl
    else baseUrl <- "http://localhost:8080/MIRACLE/readoutExport/"
    if(query$plateSecurityTokens != "") plateTokens <- str_split(query$plateSecurityTokens, "\\|")[[1]]
    else return(NULL)

    rppa.batch.load.readouts(plateSecurityTokens=plateTokens)
  })
  
  files <- reactive({
    if(is.null(input$files)) {
      # User has not uploaded a file yet
      return(NULL)
    }  
    inFile <- input$files
    return(inFile$datapath)
  })
  
  loadFiles <- function(files){
    listOfSlides <- list()
    
    for(file in files){
      currentSlide <- rppa.loadFromFile(file)
      listOfSlides[[attr(currentSlide, "slideIndex")]] <- currentSlide
    }  
    return(listOfSlides)    
  }
  
  slides <- reactive({
    if(is.null(loadAllSlides())){
      if(is.null(files())) return(NULL)
      else{
        return(loadFiles(files()))
      }
    }
    else{
      return(loadAllSlides())
    }
  })
  
  slideTitles <- function(){
    all.slides <- slides()
    if(is.null(all.slides)) return(NULL)
    titles <- c()
    slideIndices <- c()
    for(slide in all.slides){
      title <- attr(slide, "antibody")
      slideIndex <- attr(slide, "slideIndex")
      titles <- append(titles, title)
      slideIndices <- append(slideIndices, slideIndex)
    }
    names(slideIndices) <- titles
    return(slideIndices)
  }
  
  ### INPUT ELEMENTS ###
  output$slidesAvailable <- renderUI({
    all.slides <- slideTitles()
    if(is.null(all.slides)) return(NULL)
    else selectInput("selected.slides", "Choose slides", all.slides, all.slides, multiple=TRUE)    
  })
  
  output$fileUpload <- reactive({
    return(is.null(slides()))
  })
  outputOptions(output, 'fileUpload', suspendWhenHidden=FALSE)
  
  output$HKslides <- renderUI({
    all.slides <- slideTitles()
    if(is.null(all.slides)) return(NULL)
    else selectInput("selected.hk.slide", "Choose slides for housekeeping normalization", all.slides[as.integer(input$selected.slides)], multiple=TRUE)    
  })
  
  output$selectHeatmapSlide <- renderUI({
    all.slides <- slideTitles()
    if(is.null(all.slides)) return(NULL)
    else selectInput("slideSelectedForHeatmap", "Choose slide for heatmap", all.slides[as.integer(input$selected.slides)])    
  })
  
  output$selectProteinConcSlide <- renderUI({
    all.slides <- slideTitles()
    if(is.null(all.slides)) return(NULL)
    else selectInput("selectedSlideForProteinConcPlot", "Choose slides for protein concentration estimate plot", all.slides[as.integer(input$selected.slides)])    
  })
  
  output$selectSignificanceSlide <- renderUI({
    all.slides <- slideTitles()
    if(is.null(all.slides)) return(NULL)
    else selectInput("selectedSlideForSignificancePlot", "Choose slides for testing significance", all.slides[as.integer(input$selected.slides)])    
  })
  
  output$posControls <- renderUI({
    spots <- slides()[[1]]
    if(is.null(spots)) return(NULL)
    selectInput("positive.control", "Choose a positive control", setdiff(unique(spots$SpotType), "Sample"))
  })

  output$sampleSelect <- renderUI({
    spots <- slides()[[1]]
    if(is.null(spots)) return(NULL)
    sampleNames <- c(NA, as.character(sort(unique(spots$SampleName))))
    selectInput("samples", "Choose samples to include", sampleNames , selected=sampleNames, multiple=TRUE)
  })

  output$referenceSelect <- renderUI({
    spots <- slides()[[1]]
    if(is.null(spots)) return(NULL)
    selectInput("reference", "Choose reference sample (negative control)", c(as.character(sort(unique(spots$SampleName)))))
  })

  slideProperties <- function(){
    spots <- slides()[[1]]
    sort(setdiff(colnames(spots), c("vshift", "hshift", "Diameter", "Flag", "FG", "BG", "Signal", "Block", "Row", "Column", "SGADesc", "SGBDesc", "SGCDesc")))
  }

  output$selectB <- renderUI({
    selectInput("select.columns.B", "Choose horizontal sub-categories", slideProperties(), multiple=T)
  })

  output$selectA <- renderUI({
   selectInput("select.columns.A", "Choose vertical sub-categories", slideProperties(), multiple=T)
  })

  output$selectFill <- renderUI({
    selectInput("select.columns.fill", "Choose color fill categories (replicates)", slideProperties(), multiple=T)
  })
  
  output$heatmapOptions <- renderUI({
    selectInput("heatmapFill", "Select a property", c(slideProperties(), "Signal", "FG", "BG") , "Signal")
  })
  
  ### PROCESSING ###
  
  quantify <- function(slide){
    selA <- input$select.columns.A
    if(is.null(selA)) selA <- NA
    
    selB <- input$select.columns.B
    if(is.null(selB)) selB <- NA
    
    selFill <- input$select.columns.fill
    if(is.null(selFill)) selFill <- NA
    
    #selSamples <- input$select.columns.sample
    selSamples <- "SampleName"
    if(is.null(selSamples)) selFill <- NA
    switch(input$method,
           "sdc" = rppa.serialDilution(slide, select.columns.A=selA, select.columns.B=selB, select.columns.fill=selFill, make.plot=F),
           "tabus" = rppa.tabus(slide, select.columns.A=selA, select.columns.B=selB, select.columns.fill=selFill),
           "hu" = rppa.nonparam(slide, select.columns.A=selA, select.columns.B=selB, select.columns.fill=selFill),
           "supercurve" = rppa.superCurve(slide, method=input$superCurve.method, model=input$superCurve.model, make.plot=F, verbose=F, select.columns.A=selA, select.columns.B=selB, select.columns.fill=selFill))
  }
  
  processedSlides <- reactive({
    
      all.slides <- slides()
      all.slides <- all.slides[as.integer(input$selected.slides)]
      input$updateButton
      
      isolate({
        withProgress(session, min=0, max=length(all.slides), {
          counter <- 0
          setProgress(message = 'Calculation in progress',
                      detail = 'This may take a while...', value=counter)
        
          processEachSlide <- function(slide, counter){
            if(input$surfaceCorrection)
            {
              setProgress(value = counter, detail=paste("Applying surface correction to slide", attr(slide, "slideIndex")))
              slide <- rppa.surface.normalization(slide, input$positive.control)
            }
            
            if(input$quantification)
            {
              setProgress(value = counter, detail=paste("Quantifying slide", attr(slide, "slideIndex")))                          
              slide <- quantify(slide) 
            }
            
            return(slide)
          }
          
          result <- foreach(slide=all.slides) %do% {
            counter <- counter + 1
            processEachSlide(slide, counter)
          }
          names(result) <- names(all.slides)
          return(result)
        })
      })
  })
  
  normalizedSlides <- reactive({ 
    all.slides <- processedSlides()
  
    if(input$proteinLoadNormalization)
    {
      if(input$normalizationMethod == "houseKeeping"){
        if(is.null(input$selected.hk.slide)) stop("You have to select at least one slide for housekeeping normalization.")
        slidesForNormalization <- subset(all.slides, names(all.slides)%in%input$selected.hk.slide) 
      }
      else slidesForNormalization <- all.slides
      
      counter <- 1
      
      withProgress(session, min=1, max=length(all.slides), {
        setProgress(message = 'Calculation in progress',
                    detail = 'This may take a while...')
        lapply(all.slides, function(slide){
          setProgress(value = counter, detail=paste("Normalizing slide", attr(slide, "slideIndex")))
          slide <- rppa.proteinConc.normalize(slide, slidesForNormalization, method=input$normalizationMethod)
          counter <- counter + 1
          return(slide)
        })
      })
    }
    else return(all.slides)
  })
  
  formattedSlides <- reactive({
    input$updateButton
    all.slides <- normalizedSlides()
    
    isolate({
      formatEachSlide <- function(slide){
        if(!is.null(input$samples)){
          slide.readout <- attr(slide, "readout")
          slide.readout.centered <- attr(slide, "readout.centered")
          slide.readout <- slide.readout[slide[,"Sample"] %in% input$samples,]
          slide.readout.centered <- slide.readout.centered[slide[,"Sample"] %in% input$samples,]
          
          slide$Sample <- as.character(slide$Sample)
          slide <- slide[slide[,"Sample"] %in% input$samples,]
          attr(slide, "readout") <- slide.readout
          attr(slide, "readout.centered") <- slide.readout.centered
          slide$Sample <- as.factor(slide$Sample)
        }
        
        if(!is.null(input$reference) && input$normalize.to.ref.sample){
          
          slideAttr <- attributes(slide)
          slide <- rppa.normalize.to.ref.sample(slide, input$reference, each.fill=T)
          mostattributes(slide) <- slideAttr
        }
        
        return(slide)
      }    
      lapply(all.slides, formatEachSlide) 
    })  
  })
  
  
  selectedSlide <- reactive({
    all.slides <- formattedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    all.slides[[input$selectedSlideForProteinConcPlot]]
  })
  
  pairwiseCorrelations <- reactive({
    all.slides <- formattedSlides()
    
    concentrations <- foreach(slide=all.slides, .combine=cbind) %do% {
      slide$concentrations  
    } 
    colnames(concentrations) <- paste(slideTitles(), names(slideTitles()))
    
    if(input$includeReadoutInCorrelation){
      readout <- attr(all.slides[[1]], "readout")
      concentrations <- cbind(concentrations, PlateReadout=readout[,"concentrations"])
    } 
    cor(concentrations, use="pairwise.complete.obs")
  })
  
  pairwiseRawCorrelations <- reactive({
    all.slides <- slides()
    
    concentrations <- foreach(slide=all.slides, .combine=cbind) %do% {
      slide$Signal  
    } 
    colnames(concentrations) <- paste(slideTitles(), names(slideTitles()))
    cor(concentrations, use="pairwise.complete.obs")
  })
  
  dunnettsTest <- reactive({
    all.slides <- formattedSlides()
    if(is.null(input$selectedSlideForSignificancePlot)) return(NULL)
    slide <- all.slides[[input$selectedSlideForSignificancePlot]]
    if(is.null(slide)) return(NULL)
    #check that we have enough replicates
    if(is.null(slide$Fill)){
      stop("You have to choose at least one column as color fill for testing significance, since this category is used to determine replicates(check 'Show sample options' first.)")
    } else if(is.null(slide$A) && is.null(slide$B)){ checkResult <- count(slide$Sample) 
    } else if(is.null(slide$A)){ checkResult <- ddply(slide[,c("Sample", "B")], .(Sample, B), summarise, freq=length(Sample)) 
    } else if(is.null(slide$B)){ checkResult <- ddply(slide[,c("Sample", "A")], .(Sample, A), summarise, freq=length(Sample)) 
    } else { checkResult <- ddply(slide[,c("Sample", "A", "B")], .(Sample, A, B), summarise, freq=length(Sample)) 
    }
    if(min(checkResult$freq) < 2) stop("Not enough replicates for all selected samples, try a different color fill category (used to determine replicates) or exclude samples with too few replicates.")
    withProgress(session, min=1, max=5, {
      setProgress(message = 'Performing Dunnett test',
                  detail = 'a few seconds away...')
      results <- rppa.dunnett(slide=slide, referenceSample=input$reference)
      return(results)
    })   
  })
  
  ### TABLES ###

  output$proteinConcTable <- renderDataTable({
    all.slides <- formattedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    all.slides[[input$selectedSlideForProteinConcPlot]]
  })
  
  output$signDiffTable <- renderDataTable({
    dunnettsTest()  
  })
  
  ### PLOTS ###
  correlationsPlot <- function(melted.correlations, title){
    p <- qplot(x=X1, y=X2, data=melted.correlations, xlab="", ylab="", fill=value, main=title)
    p <- p + geom_tile(aes(fill = melted.correlations$value, line = 0))
    p <- p + geom_text(aes(fill = melted.correlations$value, label = round(melted.correlations$value, 2)), colour="white")
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   panel.margin = unit(0.1, "lines"), panel.margin = unit(0, "lines"), plot.margin = unit(c(1, 1, 0.5, 0.5), "lines"),
                   plot.title = element_text(size = 18), strip.background = element_rect(fill = "grey90", colour = "grey50"))
    return(p)
  }
  
  output$correlationPlot <- renderPlot({
    if(is.null(slides())) stop("No slides loaded")
    correlations <- pairwiseCorrelations()
    melted.correlations <- melt(correlations)
    print(correlationsPlot(melted.correlations, "Pearson correlation of protein concentration estimates"))
  })
  
  output$rawCorrelationPlot <- renderPlot({
    if(is.null(slides())) stop("No slides loaded")
    correlations <- pairwiseRawCorrelations()
    melted.correlations <- melt(correlations)
    print(correlationsPlot(melted.correlations, "Pearson correlation of raw signal"))
  })
  
  output$quantificationFitPlot <- renderPlot({
    if(is.null(slides())) stop("No slides loaded")
    all.slides <- formattedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    slide <- all.slides[[input$selectedSlideForProteinConcPlot]]    
    
    if(input$method == "sdc"){
      
      fittedData <- attr(slide, "fittedData")
      pairedData <- attr(slide, "pairedData")
      D <- attr(slide, "estimatedDilutionFactor")
      
      plotA <- ggplot(fittedData, aes(x=y, y=log2(x))) + labs(title="Signal vs concentration estimate plot") + ylab("Signal") + xlab("Concentration estimate") + geom_point() + geom_smooth(aes(group=1), method="loess")
      plotB <- ggplot(pairedData, aes(x=x, y=y)) + labs(title=paste("Serial Dilution Curve Fit (estimated dilution factor ", D, ")")) + xlab("Signal at next dilution step") + ylab("Signal") + geom_point() + geom_line(data=fittedData, color="blue") + geom_abline(intercept=0, slope=1, color="red")
      
      print(multiplot(plotA, plotB, cols=1))
    }
    else if(input$method == "supercurve"){
      par(mfrow(c(2,1)))
      new.fit <- attr(slide, "fit")
      plot(new.fit)
      image(new.fit)
      par(mfrow(c(1,1)))
    }
    
  })
  
  useReadoutAsFill <- function(result){  
    if(!is.null(result$B)){ 
      result$B <- paste(result$B, result$Fill, sep="|")
    } else result$B <- result$Fill
    result$B <- as.factor(result$B)
    readout <- result[,intersect(colnames(result), c("Sample", "A", "B"))]
    readout.data <- attr(result, "readout.centered")
    result <- cbind(readout, (result[,c("concentrations", "upper", "lower")] / mean(result$concentrations, na.rm=T)))
    readout <- cbind(readout, readout.data)
    result$Fill <- "Protein Conc. Est."
    readout$Fill <- "Plate Readout"

    result <- rbind(result, readout)
    if(input$normalize.to.ref.sample) result <- rppa.normalize.to.ref.sample(result, sampleReference=input$reference, each.fill=T)
      
    return(result)
  }
  
  output$proteinConcPlot <- renderPlot({
    if(is.null(slides())) stop("No slides loaded")
    all.slides <- formattedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    slide <- all.slides[[input$selectedSlideForProteinConcPlot]]
    
    if(input$compareToReadoutData) slide <- useReadoutAsFill(slide)
    
    rppa.proteinConc.plot(slide, title=attr(slide, "title"), swap=input$swap, 
                          horizontal.line=input$horizontal.line, error.bars=input$error.bars, scales=input$scales)
  })
  
  output$proteinConcOverviewPlot <- renderPlot({
    if(is.null(slides())) stop("No slides loaded")
    all.slides <- formattedSlides()
    data.protein.conc <- ldply(all.slides)
    data.protein.conc$Slide <- apply(data.protein.conc, 1, function(x){
      paste(x[1], ":", names(slideTitles()[slideTitles()==x[1]]), sep="")
    })
    
    if(input$includeReadoutInComparison){
      anySlide <- all.slides[[1]]
      readout <- cbind(anySlide[,intersect(colnames(anySlide),c("Sample", "A", "B", "Fill"))], 
                     attr(anySlide, "readout.centered"))
      if(input$normalize.to.ref.sample)
      {
        readout <- rppa.normalize.to.ref.sample(readout, sampleReference=input$reference, each.fill=T)
        readout <- readout[, setdiff(colnames(readout), c(".id", "reference"))]
      }

      readout$Slide <- "Plate Readout"
      data.protein.conc <- rbind(data.protein.conc[,-1], readout)      
    }
    
    data.protein.conc.copy <- ddply(data.protein.conc, intersect(colnames(data.protein.conc), c("Sample", "Slide", "A", "B")), summarise,
                                    concentrations = mean(concentrations, na.rm=T),
                                    upper = max(upper, na.rm=T),
                                    lower = min(lower, na.rm=T)) 

    rppa.proteinConc.plot(data.protein.conc.copy, "Protein Concentration Estimate Comparison", input$swap, input$horizontal.line, 
                          input$fill.legend, input$error.bars, input$scales, slideAsFill=T)
    
  })
  
  #create the heatmap plot
  output$heatmapPlot <- renderPlot({  
    if(is.null(slides())) stop("No slides loaded")
    if(is.null(input$slideSelectedForHeatmap)) return(NULL)
    
    rppa.plot.heatmap(slides()[[input$slideSelectedForHeatmap]], log=input$heatmapLog, fill=input$heatmapFill, discreteColorA=input$discreteColorA,
                      discreteColorB=input$discreteColorB, plotNA=input$heatmapPlotNA, palette=input$heatmapPalette)
  })
  
  output$dunnettsPlot <- renderPlot({
    testResult <- dunnettsTest()
    if(is.null(testResult)) return(NULL)
    else rppa.plot.dunnett(testResult, set.neg.control.to.one=input$sign.neg.ctrl.to.one)
  })

  ### DATA OUTPUT ###
  output$downloadData <- downloadHandler(
    filename = function() { paste(title(), input$method, '.csv', sep=' ')},
    content = function(file) {
      if(input$csvType == "CSV") write.csv(quantified(), file)
      else write.csv2(quantified(), file)
   }
  )
  
  output$downloadSignDiffData <- downloadHandler(    
    filename = function() {
      slide <- dunnettsTest()
      file <- paste(attr(slide, "title"), input$method, "_significance", sep=' ')
      switch(input$tableFileType,
             "CSV" = paste(file, ".csv", sep=""),
             "CSV2" = paste(file, ".csv", sep=""),
             "TAB" = paste(file, ".txt", sep=""),
             "XLSX" = paste(file, ".xlsx", sep=""))
    },
    content = function(file) {  
      slide <- dunnettsTest()
      switch(input$tableFileType,
             "CSV" = write.csv(slide, file),
             "CSV2" = write.csv2(slide, file),
             "TAB" = write.table(slide, file, sep="\t", row.names=F, col.names=T, quote=F),
             "XLSX" = write.xlsx(slide, file))
    }
  )
  
  output$downloadProteinConcData <- downloadHandler(    
    filename = function() {
      slide <- selectedSlide()
      file <- paste(attr(slide, "title"), input$method, sep=' ')
      switch(input$tableFileType,
             "CSV" = paste(file, ".csv", sep=""),
             "CSV2" = paste(file, ".csv", sep=""),
             "TAB" = paste(file, ".txt", sep=""),
             "XLSX" = paste(file, ".xlsx", sep=""))
    },
    content = function(file) {  
      slide <- selectedSlide()
      switch(input$tableFileType,
             "CSV" = write.csv(slide, file),
             "CSV2" = write.csv2(slide, file),
             "TAB" = write.table(slide, file, sep="\t", row.names=F, col.names=T, quote=F),
             "XLSX" = write.xlsx(slide, file))
    }
  )
}) 
