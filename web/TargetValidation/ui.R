library(shiny)
library(plotly)
library(d3heatmap)

# Define UI for random distribution application 
shinyUI(fluidPage(
    
  # Application title
  titlePanel("Target Validation"),
  
  # Sidebar with controls to select the random distribution type
  # and number of observations to generate. Note the use of the
  # br() element to introduce extra vertical spacing
  sidebarLayout(
    sidebarPanel(
      radioButtons("dataset", "Data Set:",
                   c("Ovarian TCGA" = "ov",
                     "Pancreatic TCGA" = "paad",
                     "Prostate TCGA" = "prad")),
      br(),
      
	tags$textarea(id="foo", rows=20, cols=40, "Default value"),

	submitButton("Run"),
br(),
downloadButton('downloadData', 'Download')
	
    ),

    
    # Show a tabset that includes a plot, summary, and table view
    # of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs", 
        tabPanel("Table", dataTableOutput(outputId="table")), 
        tabPanel("Waterfall", plotOutput("summary")) 
      )
    )
  )
))
