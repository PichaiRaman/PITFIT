source("../../code/aim1.R");
library(shiny)


# Define server logic for random distribution application
shinyServer(function(input, output) {
  
 
 
  # Generate a plot of the data. Also uses the inputs to build
  # the plot label. Note that the dependencies on both the inputs
  # and the data reactive expression are both tracked, and
  # all expressions are called in the sequence implied by the
  # dependency graph

  output$table <- renderTable({
    	myRes <- SigLimmaTrain(input$dataset,input$gene, thresh=.20, pvalThresh=.25, logFCThresh=1);
    	myRes[[2]];
  })

  output$plot <- renderPlot({ 
  	myRes <- SigLimmaTrain(input$dataset,input$gene, thresh=.20, pvalThresh=.25, logFCThresh=1);
	plotVolcanoTrain(myRes[[1]]);
  })
  
  # Generate a summary of the data
  output$summary <- renderPlot({
    	myRes <- SigLimmaTrain(input$dataset,input$gene, thresh=.20, pvalThresh=.25, logFCThresh=1);
	plotHeatmap(myRes[[2]], myRes[[3]]);
  })
  
  
})

