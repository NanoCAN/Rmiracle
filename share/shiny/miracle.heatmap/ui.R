library(shiny)
library(Rmiracle)
library(rCharts)

shinyUI(pageWithSidebar(
  headerPanel("Heatmap"),
  sidebarPanel(
    selectInput("log", "log-transformation?", 
                choices = c("none", "log2", "log10")),
    selectInput("fill", "Select a property", 
                choices = c("Signal", "FG", "BG","Deposition", "CellLine", "LysisBuffer", "DilutionFactor", "Inducer", "SpotType", "SpotClass", "SampleName", "SampleType", "TargetGene")),
    checkboxInput("plotNA", "Mark excluded spots", value=TRUE),
    selectInput("discreteColorA", "Select color A", 
                choices = c("darkblue", "red", "blue", "steelblue", "magenta", "yellow", "white", "green")),
    selectInput("discreteColorB", "Select color B", c("red", "darkblue", "blue", "steelblue", "magenta", "yellow", "white", "green"))
  ),
  mainPanel(
      showOutput("dynamicHeatmap", "polycharts")   
  )
))

