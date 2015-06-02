source("../code/aim1.R");
library(shiny)



# Define server logic for random distribution application
shinyServer(function(input, output) {
  
  myRes <- SigLimmaTrain(input$dataset, input$gene, thresh=.20, pvalThresh=.25, logFCThresh=1);
  
  # Generate a plot of the data. Also uses the inputs to build
  # the plot label. Note that the dependencies on both the inputs
  # and the data reactive expression are both tracked, and
  # all expressions are called in the sequence implied by the
  # dependency graph
  output$plot <- renderPlot({ 
  plotVolcanoTrain(myRes[[1]]);
  })
  
  # Generate a summary of the data
  output$summary <- renderPlot({
    plotVolcanoTrain(myRes[[1]]);
  })
  
  # Generate an HTML table view of the data
  output$table <- renderTable({
    myRes[[2]]
  })
  
})

