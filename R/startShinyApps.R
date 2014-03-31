rppa.startShinyAnalysis <- function(launch.browser=T){
	message("Starting shiny analysis app in browser")
	shiny::runApp(
	  system.file('shiny.batch.analysis',                                                    
	              package='Rmiracle'))

}
