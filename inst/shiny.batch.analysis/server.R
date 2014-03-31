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
    if(length(query$slideSecurityTokens)> 0) slideTokens <- str_split(query$slideSecurityTokens, "\\|")[[1]]
    else return(NULL)
    
    withProgress(session, min=1, max=(length(slideTokens)+1), expr={
      slides <- list()
      count <- 1
      for(token in slideTokens)
      {
        setProgress(message = 'Fetching data from MIRACLE...',
                  detail = 'Please be patient!',
                  value=count)
        currentSlide <- getSlide(baseUrl, token)
        slides[[attr(currentSlide, "slideIndex")]] <- currentSlide
        count <- count + 1
      }
      return(slides)
    })
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
    else selectInput("selected.hk.slide", "Choose slides for housekeeping normalization", all.slides, all.slides[[1]], multiple=TRUE)    
  })
  
  output$selectHeatmapSlide <- renderUI({
    all.slides <- slideTitles()
    if(is.null(all.slides)) return(NULL)
    else selectInput("slideSelectedForHeatmap", "Choose slide for heatmap", all.slides, all.slides[[1]])    
  })
  
  output$selectProteinConcSlide <- renderUI({
    all.slides <- slideTitles()
    if(is.null(all.slides)) return(NULL)
    else selectInput("selectedSlideForProteinConcPlot", "Choose slides for protein concentration estimate plot", all.slides, all.slides[[1]])    
  })
  
  output$selectSignificanceSlide <- renderUI({
    all.slides <- slideTitles()
    if(is.null(all.slides)) return(NULL)
    else selectInput("selectedSlideForSignificancePlot", "Choose slides for testing significance", all.slides, all.slides[[1]])    
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

  output$selectA <- renderUI({
    selectInput("select.columns.A", "Choose horizontal sub-categories", slideProperties(), multiple=T)
  })

  output$selectB <- renderUI({
   selectInput("select.columns.B", "Choose vertical sub-categories", slideProperties(), multiple=T)
  })

  output$selectFill <- renderUI({
    selectInput("select.columns.fill", "Choose color fill categories (replicates)", slideProperties(), multiple=T)
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
           "hu" = rppa.nonparam(slide, , select.columns.A=selA, select.columns.B=selB, select.columns.fill=selFill),
           "supercurve" = rppa.superCurve(slide, , select.columns.A=selA, select.columns.B=selB, select.columns.fill=selFill))
  }
  
  processedSlides <- reactive({
    all.slides <- slides()
    counter <- 1
    
    withProgress(session, min=1, max=length(all.slides)*2, {
      setProgress(message = 'Calculation in progress',
                  detail = 'This may take a while...')
      
      processEachSlide <- function(slide){
        if(input$surfaceCorrection)
        {
          setProgress(value = counter, detail=paste("Applying surface correction to slide", attr(slide, "slideIndex")))
          slide <- rppa.surface.normalization(slide, input$positive.control)
        }
        counter <- counter + 1
        if(input$quantification)
        {
          setProgress(value = counter, detail=paste("Quantifying slide", attr(slide, "slideIndex")))                          
          slide <- quantify(slide) 
        }
        return(slide)
      }
      
      lapply(all.slides, processEachSlide) 
    })
  })

  normalizedSlides <- reactive({    
    all.slides <- processedSlides()
    if(input$proteinLoadNormalization)
    {
      if(input$normalizationMethod == "houseKeeping"){
        if(is.null(input$selected.hk.slide)) return(NULL)
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
  
  selectedSlide <- reactive({
    all.slides <- normalizedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    all.slides[[input$selectedSlideForProteinConcPlot]]
  })
  
  pairwiseCorrelations <- reactive({
    all.slides <- normalizedSlides()
    concentrations <- foreach(slide=all.slides, .combine=cbind) %do% {
      slide$concentrations  
    } 
    colnames(concentrations) <- paste(slideTitles(), names(slideTitles()))
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
    all.slides <- normalizedSlides()
    if(is.null(input$selectedSlideForSignificancePlot)) return(NULL)
    slide <- all.slides[[input$selectedSlideForSignificancePlot]]
    if(is.null(slide)) return(NULL)
    #check that we have enough replicates
    if(is.null(slide$Fill)){
      stop("You have to choose at least one column as color fill for testing significance, since this category is used to determine replicates(check 'Show sample options' first.)")
    } else if(is.null(slide$A) && is.null(slide$B)){ checkResult <- ddply(slide, .(Sample), summarise, freq=length(Sample)) 
    } else if(is.null(slide$A)){ checkResult <- ddply(slide, .(Sample, B), summarise, freq=length(Sample)) 
    } else if(is.null(slide$B)){ checkResult <- ddply(slide, .(Sample, A), summarise, freq=length(Sample)) 
    } else { checkResult <- ddply(slide[,c("Sample", "A", "B")], .(Sample, A, B), summarise, freq=length(Sample)) 
    }
    if(min(checkResult$freq) < 2) stop("Not enough replicates for all selected samples, try a different color fill category (used to determine replicates) or exclude samples with too few replicates.")
    withProgress(session, min=1, max=5, {
      setProgress(message = 'Performing Dunnett test',
                  detail = 'a few seconds away...')
      results <- rppa.dunnett(slide=slide, referenceSample=input$reference, sample.subset=input$samples)
      return(results)
    })   
  })
  
  ### TABLES ###

  output$proteinConcTable <- renderDataTable({
    all.slides <- normalizedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    all.slides[[input$selectedSlideForProteinConcPlot]]
  })
  
  ### PLOTS ###
  correlationsPlot <- function(melted.correlations, title){
    p <- qplot(x=X1, y=X2, data=melted.correlations, xlab="", ylab="", fill=value, main=title)
    p <- p + geom_tile(aes(fill = melted.correlations$value, line = 0))
    p <- p + geom_text(aes(fill = melted.correlations$value, label = round(melted.correlations$value, 2)))
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
    all.slides <- normalizedSlides()
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
    
  })
  
  output$proteinConcPlot <- renderPlot({
    if(is.null(slides())) stop("No slides loaded")
    all.slides <- normalizedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    slide <- all.slides[[input$selectedSlideForProteinConcPlot]]
    reference <- input$reference
    if(is.null(reference)) reference <- NA
    
    samples <- input$samples
    if(is.null(samples)) samples <- levels(slide$Sample)
    
    rppa.proteinConc.plot(slide, title=attr(slide, "title"), swap=input$swap, 
                          horizontal.line=input$horizontal.line, error.bars=input$error.bars, scales=input$scales,
                          sample.subset=samples, reference=reference)
  })
  
  output$proteinConcOverviewPlot <- renderPlot({
    if(is.null(slides())) stop("No slides loaded")
    all.slides <- normalizedSlides()
    #slides <- lapply(slides, function(x){ x$Slide <- attr(x, "antibody"); return(x)})
    data.protein.conc <- ldply(all.slides)
    data.protein.conc$Slide <- apply(data.protein.conc, 1, function(x){
      paste(x[1], ":", names(slideTitles()[slideTitles()==x[1]]), sep="")
    })

    reference <- input$reference
    if(is.null(reference)) reference <- NA
    
    samples <- input$samples
    if(is.null(samples)) samples <- levels(slide$Sample)
    
    data.protein.conc.copy <- data.protein.conc
    
    normalizeFill <- function(data.protein.conc){
      ddply(data.protein.conc, intersect(colnames(data.protein.conc), c("Fill", "A", "B")), transform, 
            concentrations = concentrations / mean(concentrations, na.rm=T),
            upper = upper / mean(concentrations, na.rm=T),
            lower = lower / mean(concentrations, na.rm=T)) 
    }
    
    data.protein.conc.copy <- ddply(normalizeFill(data.protein.conc), intersect(colnames(data.protein.conc), c("Sample", "Slide", "A", "B")), summarise,
                                    concentrations = mean(concentrations, na.rm=T),
                                    upper = max(upper, na.rm=T),
                                    lower = min(lower, na.rm=T)) 
    each.A <- F
    each.B <- F
    specific.A.copy <- NULL
    specific.B.copy <- NULL
    each.fill <- T
    fill.legend <- T
    
    #if(duplicate.na)
    #{
      #data.protein.conc.copy <- rppa.duplicate.nas(data.protein.conc.copy)
    #}    
    #rppa.proteinConc.plot(data.protein.conc.copy, slideAsFill=T)
    rppa.proteinConc.plot(data.protein.conc.copy, "Protein Concentration Estimate Comparison", input$swap, input$horizontal.line, 
                          input$fill.legend, input$error.bars, input$scales, samples, reference, slideAsFill=T)
    
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