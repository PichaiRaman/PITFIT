source("../../code/aim-2/aim-2-worker.R");
library(shiny)
library(d3heatmap)
library(plotly);

globalResult <- "";
globalVolcano <- "";
# Define server logic for random distribution application
shinyServer(function(input, output, session) {
  
 
  updateSelectizeInput(session, 'gene', choices = c("FOXM1", "MYC", "E2F1"), server = TRUE)
 
  # Generate a plot of the data. Also uses the inputs to build
  # the plot label. Note that the dependencies on both the inputs
  # and the data reactive expression are both tracked, and
  # all expressions are called in the sequence implied by the
  # dependency graph

  output$table <- renderDataTable({
  
      input$submit # Re-run when button is clicked

    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))

    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Running Analysis", value = 0)

    # Number of times we'll go through the loop
    n <- 50

    for (i in 1:n) {
      # Each time through the loop, add another row of data. This is
      # a stand-in for a long-running computation.
      dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))

      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", i))

      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
    }

  	myRes <- pitfitAnalyzeAim2(input$dataset,input$gene);
globalResult <<- myRes; 
myRes;
  })

 
  output$downloadData <- downloadHandler(
    filename = function() { paste('result', '.csv', sep='') },
    content = function(file) {
      write.csv(globalResult, file, row.names=F)
    })


  
})
