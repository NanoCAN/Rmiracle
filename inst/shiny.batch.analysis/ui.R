library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("Rmiracle"),
  sidebarPanel(
    tags$head( 
      tags$link(rel="stylesheet", type="text/css", 
                href="rmiracle.css") 
    ),
    conditionalPanel("output.fileUpload",
      fileInput("files", "File data", multiple=TRUE)
    ),    
    actionButton("updateButton", "Update"),
    br(),br(),
    uiOutput("slidesAvailable"),
    checkboxInput("surfaceCorrection", "Use positive control spots correcting staining bias?", FALSE),
    conditionalPanel(condition= "input.surfaceCorrection",
      uiOutput("posControls")
    ),
    checkboxInput("quantification", "Quantification of dilution series?", TRUE),
    conditionalPanel(condition = "input.quantification",   
      selectInput("method", "Quantification Method:",
                choices = c("Serial Dilution Curve"="sdc", "SuperCurve"="supercurve", "Tabus"="tabus", "Non-parametric"="hu"))),
      conditionalPanel(condition="input.method=='supercurve'",
                strong("Warning: Computation with this method can be time intensive or fail, probably due to convergence issues with few dilution steps."),
      selectInput("superCurve.method", "Select a fit method", c("nls", "nlrob", "nlrq"), "nlrq"),
      selectInput("superCurve.model", "Select a response curve model", c("logistic", "loess", "cobs"), "cobs")),
      conditionalPanel(condition="input.method=='hu'", strong("Warning: Computation with this method is time intensive")),
    checkboxInput("estimateNormalization", "Apply quantile normalization?", FALSE),
    conditionalPanel(condition= "input.estimateNormalization",
                     selectInput("estimateNormMethod", "Normalize:", choices = c("Affy"="affyQuant", "Aroma light"="aromaQuant"))                  
    ),
    checkboxInput("calcNormFactors", "Calculate normalization factors for each slide?", FALSE),
    conditionalPanel(condition= "input.calcNormFactors",
                     selectInput("calcNormMethod", "Normalize:", choices = c("Trimmed Mean"="TMM", "Relative log expression"="RLE", "upperquartile"="upper quartile"))                     
    ),
    checkboxInput("proteinLoadNormalization", "Normalize for total protein amount?", FALSE),
    conditionalPanel(condition= "input.proteinLoadNormalization",
                     selectInput("normalizationMethod", "Normalization Method:", choices = c("House keeping proteins"="houseKeeping", "Median Loading"="medianLoading", "Variable Slope"="variableSlope")),
                     conditionalPanel(condition="input.normalizationMethod=='houseKeeping'",
                        uiOutput("HKslides"))
    ),
    #selectInput("normalizeDepositions", "Merge depositions?", choices= c("No"="None", "Mean"="mean", "Linear Regression"="LinReg"), "None"),
    checkboxInput("groupingOptions", "Show sample options?", FALSE),
    conditionalPanel(condition = "input.groupingOptions",
                     uiOutput("sampleSelect"),
                     checkboxInput("normalize.to.ref.sample", "Normalize to reference sample?", value=FALSE),              
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
      selectInput("tableFileType", "Select file type for table downloads", choices=c("CSV", "CSV2", "TAB"), selected="TAB")
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Heatmap", 
        uiOutput("selectHeatmapSlide"),
        checkboxInput("heatmapPlotOptions", "Plot options", FALSE),
        conditionalPanel(condition = "input.heatmapPlotOptions", 
                         fluidRow(
                           column(3,
                                  selectInput("heatmapLog", "log-transformation?", 
                                              choices = c("none", "log2", "log10")),
                                  uiOutput("heatmapOptions"),
                                  checkboxInput("heatmapPlotNA", "Mark excluded spots", value=TRUE)
                           ), 
                           column(3, offset = 1,
                                  selectInput("heatmapPalette", "Select color palette for categorical variables", choices=c("Set1", "Set2", "Set3", 
                                                "Accent", "Dark2", "Paired", "Pastel1", "Pastel2"), selected="Set1"),
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
               checkboxInput("compareToReadoutData", "Compare to readout data?", value=FALSE),
               plotOutput("proteinConcPlot"),
               dataTableOutput("proteinConcTable"), 
               downloadButton('downloadProteinConcData', 'Download'),
               plotOutput("quantificationFitPlot")),
      tabPanel("Comparison", 
          checkboxInput("includeReadoutInComparison", "Include readout data?", value=FALSE),     
          plotOutput("proteinConcOverviewPlot")),
      tabPanel("Correlation",     
          checkboxInput("includeReadoutInCorrelation", "Include readout in correlation?", value=FALSE),                    
          plotOutput("correlationPlot"),
          plotOutput("rawCorrelationPlot")),
      tabPanel("Significance",
          checkboxInput("sign.neg.ctrl.to.one", "Set reference (negative control) to 100%?", value=F),
          uiOutput("selectSignificanceSlide"),                         
          plotOutput("dunnettsPlot"),
          dataTableOutput("signDiffTable"), 
          downloadButton('downloadSignDiffData', 'Download'))
    )
  )
))
