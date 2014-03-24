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
    slideTokens <- str_split(query$slideSecurityTokens, "\\|")[[1]]
    
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
  
  slideTitles <- function(){
    slides <- loadAllSlides()
    titles <- c()
    slideIndices <- c()
    for(slide in slides){
      title <- attr(slide, "antibody")
      slideIndex <- attr(slide, "slideIndex")
      titles <- append(titles, title)
      slideIndices <- append(slideIndices, slideIndex)
    }
    names(slideIndices) <- titles
    return(slideIndices)
  }
  
  ### INPUT ELEMENTS ###
  output$slides <- renderUI({
    slides <- slideTitles()
    selectInput("selected.slides", "Choose slides", slides, slides, multiple=TRUE)    
  })
  
  output$HKslides <- renderUI({
    slides <- slideTitles()
    selectInput("selected.hk.slide", "Choose slides for housekeeping normalization", slides, slides[[1]], multiple=TRUE)    
  })
  
  output$selectHeatmapSlide <- renderUI({
    slides <- slideTitles()
    selectInput("slideSelectedForHeatmap", "Choose slide for heatmap", slides, slides[[1]])    
  })
  
  output$selectProteinConcSlide <- renderUI({
    slides <- slideTitles()
    selectInput("selectedSlideForProteinConcPlot", "Choose slides for protein concentration estimate plot", slides, slides[[1]])    
  })
  
  output$posControls <- renderUI({
    spots <- loadAllSlides()[[1]]
    selectInput("positive.control", "Choose a positive control", setdiff(unique(spots$SpotType), "Sample"))
  })

  output$sampleSelect <- renderUI({
    spots <- loadAllSlides()[[1]]
    sampleNames <- c(NA, as.character(sort(unique(spots$SampleName))))
    selectInput("samples", "Choose samples to include", sampleNames , selected=sampleNames, multiple=TRUE)
  })

  output$referenceSelect <- renderUI({
    spots <- loadAllSlides()[[1]]
    selectInput("reference", "Choose reference sample", c(as.character(sort(unique(spots$SampleName)))))
  })

  slideProperties <- function(){
    spots <- loadAllSlides()[[1]]
    sort(setdiff(colnames(spots), c("vshift", "hshift", "Diameter", "Flag", "PlateLayout", "PlateRow", "PlateCol", "FG", "BG", "Signal", "Block", "Row", "Column")))
  }

  output$selectA <- renderUI({
    selectInput("select.columns.A", "Choose horizontal sub-categories", slideProperties(), multiple=T)
  })

  output$selectB <- renderUI({
   selectInput("select.columns.B", "Choose vertical sub-categories", slideProperties(), multiple=T)
  })

  output$selectFill <- renderUI({
    selectInput("select.columns.fill", "Choose color fill category", slideProperties())
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
    slides <- loadAllSlides()
    counter <- 1
    
    withProgress(session, min=1, max=length(slides)*2, {
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
      
      lapply(slides, processEachSlide) 
    })
  })

  normalizedSlides <- reactive({    
    slides <- processedSlides()
    if(input$proteinLoadNormalization)
    {
      if(input$normalizationMethod == "houseKeeping"){
        if(is.null(input$selected.hk.slide)) return(NULL)
        slidesForNormalization <- subset(slides, names(slides)%in%input$selected.hk.slide)
      }
      else slidesForNormalization <- slides
      
      counter <- 1
      
      withProgress(session, min=1, max=length(slides), {
        setProgress(message = 'Calculation in progress',
                    detail = 'This may take a while...')
        lapply(slides, function(slide){
          setProgress(value = counter, detail=paste("Normalizing slide", attr(slide, "slideIndex")))
          slide <- rppa.proteinConc.normalize(slide, slidesForNormalization, method=input$normalizationMethod)
          counter <- counter + 1
          return(slide)
        })
      })
    }
    else return(slides)
  })
  
  pairwiseCorrelations <- reactive({
    slides <- normalizedSlides()
    concentrations <- foreach(slide=slides, .combine=cbind) %do% {
      slide$concentrations  
    } 
    colnames(concentrations) <- paste(slideTitles(), names(slideTitles()))
    cor(concentrations, use="pairwise.complete.obs")
  })
  
  pairwiseRawCorrelations <- reactive({
    slides <- loadAllSlides()
    concentrations <- foreach(slide=slides, .combine=cbind) %do% {
      slide$Signal  
    } 
    colnames(concentrations) <- paste(slideTitles(), names(slideTitles()))
    cor(concentrations, use="pairwise.complete.obs")
  })
  
  ### TABLES ###

  output$proteinConcTable <- renderDataTable({
    slides <- normalizedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    slides[[input$selectedSlideForProteinConcPlot]]
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
    correlations <- pairwiseCorrelations()
    melted.correlations <- melt(correlations)
    print(correlationsPlot(melted.correlations, "Pearson correlation of protein concentration estimates"))
  })
  
  output$rawCorrelationPlot <- renderPlot({
    correlations <- pairwiseRawCorrelations()
    melted.correlations <- melt(correlations)
    print(correlationsPlot(melted.correlations, "Pearson correlation of raw signal"))
  })
  
  output$quantificationFitPlot <- renderPlot({
    slides <- normalizedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    slide <- slides[[input$selectedSlideForProteinConcPlot]]    
    
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
    slides <- normalizedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    slide <- slides[[input$selectedSlideForProteinConcPlot]]
    reference <- input$reference
    if(is.null(reference)) reference <- NA
    
    samples <- input$samples
    if(is.null(samples)) return(NULL)
    rppa.proteinConc.plot(slide, title=attr(slide, "title"), swap=input$swap, 
                          horizontal.line=input$horizontal.line, error.bars=input$error.bars, scales=input$scales,
                          sample.subset=samples, reference=reference)
    #rppa.proteinConc.plot(slide, title=attr(slide, "title"), swap=input$swap, 
    #  horizontal.line=input$horizontal.line, error.bars=input$error.bars, scales=input$scales, sample.subset=samples, reference=reference,
    #  each.A=T, each.B=T)
  })
  
  output$proteinConcOverviewPlot <- renderPlot({
    slides <- normalizedSlides()
    #slides <- lapply(slides, function(x){ x$Slide <- attr(x, "antibody"); return(x)})
    data.protein.conc <- ldply(slides)
    data.protein.conc$Slide <- apply(data.protein.conc, 1, function(x){
      paste(x[1], ":", names(slideTitles()[slideTitles()==x[1]]), sep="")
    })

    reference <- input$reference
    if(is.null(reference)) reference <- NA
    
    samples <- input$samples
    if(is.null(samples)) return(NULL)
    
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
    if(is.null(input$slideSelectedForHeatmap)) return(NULL)
    rppa.plot.heatmap(loadAllSlides()[[input$slideSelectedForHeatmap]], log=input$heatmapLog, fill=input$heatmapFill, discreteColorA=input$discreteColorA,
                      discreteColorB=input$discreteColorB, plotNA=input$heatmapPlotNA, palette=input$heatmapPalette)
  })

  ### DATA OUTPUT ###
  output$downloadData <- downloadHandler(
    filename = function() { paste(title(), input$method, '.csv', sep=' ')},
    content = function(file) {
      if(input$csvType == "CSV") write.csv(quantified(), file)
      else write.csv2(quantified(), file)
   }
  )
  
  selectedSlide <- reactive({
    slides <- normalizedSlides()
    if(is.null(input$selectedSlideForProteinConcPlot)) return(NULL)
    slides[[input$selectedSlideForProteinConcPlot]]
  })
  
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
