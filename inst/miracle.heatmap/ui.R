library(shiny)
library(Rmiracle)
#library(rCharts)

shinyUI(fluidPage(
  
    title = "Heatmap",
    plotOutput("heatmapPlot", height="500px"),
    hr(),
    # Add custom CSS & Javascript;
    tagList(
      tags$head(
        tags$link(rel="stylesheet", type="text/css",href="style.css"),
        tags$script(type="text/javascript", src = "busy.js")
      )
    ),
    div(class = "busy",  
        p("Calculation in progress"), 
        img(src="http://imageshack.us/a/img827/4092/ajaxloaderq.gif")
    ),
    fluidRow(
      column(3,
        selectInput("log", "log-transformation?", 
                choices = c("none", "log2", "log10")),
        selectInput("fill", "Select a property", 
                choices = c("Signal", "FG", "BG","Deposition", "CellLine", "LysisBuffer", "DilutionFactor", "Inducer", "SpotType", "SpotClass", "SampleName", "SampleType", "TargetGene")),
        checkboxInput("plotNA", "Mark excluded spots", value=TRUE)
      ), 
      column(3, offset = 1,
        selectInput("discreteColorA", "Select color A", 
                choices = c("darkblue", "red", "blue", "steelblue", "magenta", "yellow", "white", "green")),
        selectInput("discreteColorB", "Select color B", c("red", "darkblue", "blue", "steelblue", "magenta", "yellow", "white", "green")),
        selectInput("heatmapPalette", "Select color palette for categorical variables", choices=c("Set1", "Set2", "Set3", 
                                                                                                  "Accent", "Dark2", "Paired", "Pastel1", "Pastel2"), selected="Set1")
      )
    )
))

