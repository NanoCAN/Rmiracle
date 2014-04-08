library(shiny)
library(Rmiracle)
#libary(rCharts)
library(plyr)

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
  
  #create the heatmap plot
  output$heatmapPlot <- renderPlot({    
    rppa.plot.heatmap(currentSlide(), log=input$log, fill=input$fill, discreteColorA=input$discreteColorA,
                      discreteColorB=input$discreteColorB, plotNA=input$plotNA, palette=input$heatmapPalette)
  })
  
#   #dynamic heatmap
#   output$dynamicHeatmap <- renderChart({
#     spots <- currentSlide()
#     
#     tooltip <- paste0("function(item){return 'Sample: ' + item.SampleName")
#     for(i in c("Block", "Row", "Column", "Deposition", "CellLine", "NumberOfCellsSeeded", "Treatment", "LysisBuffer", "Signal", "FG", "BG")){
#       x <- paste0(" + '\n' + '", i, ": ' + item.", i)
#       tooltip <- paste0(tooltip, x)
#     }
#     tooltip <- paste0(tooltip, "}")  
#     
#     spots <- spots[with(spots, order(Block, Column, Row)),] #assure correct order
#     p <- rPlot(x="bin(Column,1)", y="bin(Row,1)", color = 'Signal', data = spots, type = 'tile', height = 800, tooltip=tooltip) 
#     p$guides(x=list(title="Column"), y=list(title="Row"))
#     #colorsForHeatmap <- paste("{color: {scale: {type: gradient, lower: ", input$discreteColorA, ", upper: ", input$discreteColorB, "}}}", sep="")
#     #cat(colorsForHeatmap)
#     #p$guides(list(color=list(scale=list(type="gradient", lower="steelblue", upper="white"))))
#     p$guides("{color: {scale: {type: gradient, lower: white, upper: steelblue}}}") 
#     
#     p$facet(var='Block', type='wrap', rows=4)
#     p$addParams(dom='dynamicHeatmap', width="1024", height="800")
#     return(p)
#   })
}) 