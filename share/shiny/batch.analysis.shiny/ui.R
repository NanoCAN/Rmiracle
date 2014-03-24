library(shiny)
library(shinyIncubator)

shinyUI(pageWithSidebar(
  headerPanel("Rmiracle"),
  sidebarPanel(    
    uiOutput("slides"),
    checkboxInput("surfaceCorrection", "Use positive control spots correcting staining bias?", FALSE),
    conditionalPanel(condition= "input.surfaceCorrection",
      uiOutput("posControls")
    ),
    checkboxInput("quantification", "Quantification of dilution series?", TRUE),
    conditionalPanel(condition = "input.quantification",   
      selectInput("method", "Quantification Method:",
                choices = c("Serial Dilution Curve"="sdc", "SuperCurve"="supercurve", "Tabus"="tabus", "Non-parametric"="hu"))),
    checkboxInput("proteinLoadNormalization", "Normalize for total protein amount?", FALSE),
    conditionalPanel(condition= "input.proteinLoadNormalization",
                     selectInput("normalizationMethod", "Normalization Method:", choices = c("House keeping proteins"="houseKeeping", "Median Loading"="medianLoading", "Variable Slope"="variableSlope")),
                     conditionalPanel(condition="input.normalizationMethod=='houseKeeping'",
                        uiOutput("HKslides"))
    ),
    checkboxInput("groupingOptions", "Show sample options?", FALSE),
    conditionalPanel(condition = "input.groupingOptions",
                     uiOutput("sampleSelect"),
                     uiOutput("referenceSelect"),
                     uiOutput("selectA"),
                     uiOutput("selectB"),
                     uiOutput("selectFill")),
    checkboxInput("plottingOptions", "Show plotting options?", FALSE),
    conditionalPanel(condition = "input.plottingOptions",
      selectInput("scales", "Should results be scaled?", choices = c("Scale each group separately"="fixed", "Scale only vertical"="free_y", "Scale only horizontal"="free_x", "No scaling"="free"), "free"),
      checkboxInput("swap", "Swap categories?", value=FALSE),
      checkboxInput("horizontal.line", "Draw a horizontal line through 1", value=FALSE),
      checkboxInput("error.bars", "Include error bars", value=TRUE)),
    #submitButton("update plot"),
      selectInput("tableFileType", "Select file type for table downloads", choices=c("CSV", "CSV2", "TAB", "XLSX"), selected="TAB")
  ),
  mainPanel(
    progressInit(),
    tabsetPanel(
      tabPanel("Heatmap", 
        uiOutput("selectHeatmapSlide"),
        checkboxInput("heatmapPlotOptions", "Plot options", FALSE),
        conditionalPanel(condition = "input.heatmapPlotOptions", 
                         fluidRow(
                           column(3,
                                  selectInput("heatmapLog", "log-transformation?", 
                                              choices = c("none", "log2", "log10")),
                                  selectInput("heatmapFill", "Select a property", 
                                              choices = c("Signal", "FG", "BG","Deposition", "CellLine", "LysisBuffer", "DilutionFactor", "Inducer", "SpotType", "SpotClass", "SampleName", "SampleType", "TargetGene")),
                                  checkboxInput("heatmapPlotNA", "Mark excluded spots", value=TRUE)
                           ), 
                           column(3, offset = 1,
                                  selectInput("heatmapPalette", "Select color palette for categorical variables", choices=c(NA, "Set1", "Set2", "Set3", 
                                                "Accent", "Dark2", "Paired", "Pastel1", "Pastel2"), selected=NA),
                                  selectInput("discreteColorA", "Select color A", 
                                              choices = c("darkblue", "red", "blue", "steelblue", "magenta", "yellow", "white", "green")),
                                  selectInput("discreteColorB", "Select color B", c("red", "darkblue", "blue", "steelblue", "magenta", "yellow", "white", "green"))
                           )
                         )
        ),
        plotOutput("heatmapPlot",height=800)
      ),
      tabPanel("Protein Concentration Estimates",
               uiOutput("selectProteinConcSlide"),          
               plotOutput("proteinConcPlot"),
               dataTableOutput("proteinConcTable"), 
               downloadButton('downloadProteinConcData', 'Download'),
               plotOutput("quantificationFitPlot")),
      tabPanel("Comparison", plotOutput("proteinConcOverviewPlot")),
      tabPanel("Correlation",     
          plotOutput("correlationPlot"),
          plotOutput("rawCorrelationPlot"))
    )
    #plotOutput("proteinConcPlot", height=800)
  )
))
