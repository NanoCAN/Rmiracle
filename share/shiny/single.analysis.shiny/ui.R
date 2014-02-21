library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("Rmiracle interactive analysis with shiny"),
  sidebarPanel(
    checkboxInput("surface.correction", "Use positive control spots for surface correction?", value=FALSE),
    uiOutput("posControls"),
    uiOutput("normSlides"),
    selectInput("normalizationMethod", "Normalization Method:", choices = c("medianLoading", "variableSlope")),
    selectInput("method", "Quantification Method:",
                choices = c("Serial Dilution Curve" , "SuperCurve", "Tabus", "Hu non-parametric")),
    selectInput("scales", "Should results be scaled?", choices = c("fixed", "free_y", "free_x", "free")),
    uiOutput("sampleSelect"),
    uiOutput("referenceSelect"),
    uiOutput("selectA"),
    uiOutput("selectB"),
    uiOutput("selectFill"),
    checkboxInput("swap", "Swap categories?", value=FALSE),
    checkboxInput("horizontal.line", "Draw a horizontal line through 1", value=FALSE),
    checkboxInput("error.bars", "Include error bars", value=TRUE),
    submitButton("update plot"),
    selectInput("csvType", "Select between comma separated (CSV) and semicolon separated (CSV2)", choices=c("CSV", "CSV2")),
    downloadButton('downloadData', 'download results')
  ),
  mainPanel(
    plotOutput("proteinConcPlot", height=800)
  )
))
