library(shiny)
library(stringr)
library(Rmiracle)

shinyServer(function(input, output, session) {
  
  # Parse the GET query string
  output$queryText <- renderText({
    query <- parseQueryString(session$clientData$url_search)
    
    # Return a string with key-value pairs
    paste(names(query), query, sep = "=", collapse=", ")
  })
  
  getSlide <- function(baseUrl, securityToken){
    rppa.load(baseUrl=baseUrl, securityToken=securityToken)
  }  
  
  currentSlide <- reactive({
    query <- parseQueryString(session$clientData$url_search)
    if(length(query$baseUrl > 0)) baseUrl <- query$baseUrl
    else baseUrl <- "http://localhost:8080/MIRACLE/spotExport/"
    getSlide(baseUrl, query$securityToken)
  })
  
  normalizationSlides <- reactive({
    query <- parseQueryString(session$clientData$url_search)
    if(length(query$baseUrl > 0)) baseUrl <- query$baseUrl
    else baseUrl <- "http://localhost:8080/MIRACLE/spotExport/"
    normalizationTokens <- str_split(query$normalizationTokens, "\\|")[[1]]
    normalizationSlides <- list()
    
    count <- 1
    for(token in normalizationTokens)
    {
      slide <- getSlide(baseUrl, token)
      normalizationSlides[[count]] <- slide
      count <- count + 1
    }
    return(normalizationSlides)
  })

  title <- function(){
    spots <-  currentSlide()
    attr(spots, "title")
  }
  
  normalizationTitles <- function(){
    normSlides <- normalizationSlides()
    normTitles <- c()
    for(slide in normSlides){
      title <- attr(slide, "title")
      normTitles <- append(normTitles, title)
    }
    return(normTitles)
  }
  
  output$normSlides <- renderUI({
    normSlides <- normalizationTitles()
    print(normSlides)
    selectInput("selected.norm.slides", "Choose slides for normalization", normSlides, multiple=TRUE)
  })
  
  output$posControls <- renderUI({
    spots <- currentSlide()
    selectInput("positive.control", "Choose a positive control", setdiff(unique(spots$SpotType), "Sample"))
  })

  output$sampleSelect <- renderUI({
    spots <- currentSlide()
    sampleNames <- c(NA, as.character(sort(unique(spots$SampleName))))
    selectInput("samples", "Choose samples to include", sampleNames , selected=sampleNames, multiple=TRUE)
  })

  output$referenceSelect <- renderUI({
    spots <- currentSlide()
    selectInput("reference", "Choose reference sample", c(NA, as.character(sort(unique(spots$SampleName)))))
  })

  slideProperties <- function(){
    spots <- currentSlide()
    sort(setdiff(colnames(currentSlide()), c("vshift", "hshift", "Diameter", "Flag", "PlateLayout", "PlateRow", "PlateCol", "FG", "BG", "Signal", "Block", "Row", "Column")))
  }

  output$selectA <- renderUI({
    spots <- currentSlide()
    slProps <- slideProperties()
    selectInput("select.columns.A", "Choose horizontal sub-category", slProps, multiple=T, selected="CellLine")
  })

  output$selectB <- renderUI({
   spots <- currentSlide()
   slProps <- slideProperties()
   selectInput("select.columns.B", "Choose vertical sub-category", slProps, multiple=T, select="Treatment")
  })

  output$selectFill <- renderUI({
    spots <- currentSlide()
    slProps <- slideProperties()
    selectInput("select.columns.fill", "Choose color fill category", slProps, select="NumberOfCellsSeeded")
  })
  
  preProcessed <- reactive({
    if(input$surface.correction) 
      rppa.surface.normalization(currentSlide(), input$positive.control)
    else currentSlide()
  })
  
  preProcessedNormalizationSlides <- reactive({
    normSlides <- normalizationSlides()
    normTitles <- normalizationTitles()
    selectedForNorm <- input$selected.norm.slides
    selectedNormSlides <- normSlides[normTitles %in% selectedForNorm]
    if(input$surface.correction){
      counter <- 1
      for(slide in selectedNormSlides){
        selectedNormSlides[[counter]] <- rppa.surface.normalization(slide, input$positive.control)
      }
    }
    selectedNormSlides
  })

  quantified <- reactive({
    switch(input$method,
           "Serial Dilution Curve" = rppa.serialDilution(preProcessed(), select.columns.A=input$select.columns.A, select.columns.B=input$select.columns.B, select.columns.fill=input$select.columns.fill),
           "Tabus" = rppa.tabus(preProcessed(), select.columns.A=input$select.columns.A, select.columns.B=input$select.columns.B, select.columns.fill=input$select.columns.fill),
           "Hu non-parametric" = rppa.nonparam(preProcessed(), select.columns.A=input$select.columns.A, select.columns.B=input$select.columns.B, select.columns.fill=input$select.columns.fill),
           "SuperCurve" = rppa.superCurve(preProcessed(), select.columns.A=input$select.columns.A, select.columns.B=input$select.columns.B, select.columns.fill=input$select.columns.fill))
  })
  
  normalized <- reactive({
    rppa.proteinConc.normalize(preProcessed(), preProcessedNormalizationSlides(), method=input$normalizationMethod)
  })
  
  output$proteinConcPlot <- renderPlot({
    reference <- input$reference
    if(length(reference) == 0) reference <- NA
    rppa.proteinConc.plot(quantified(), title=title(), swap=input$swap, 
      horizontal.line=input$horizontal.line, error.bars=input$error.bars, scales=input$scales, sample.subset=input$samples, reference=reference,
      each.A=T, each.B=T)
  })

  output$downloadData <- downloadHandler(
    filename = function() { paste(title(), input$method, '.csv', sep=' ')},
    content = function(file) {
      if(input$csvType == "CSV") write.csv(quantified(), file)
      else write.csv2(quantified(), file)
   }
  )
}) 
