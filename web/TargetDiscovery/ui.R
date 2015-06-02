library(shiny)

# Define UI for random distribution application 
shinyUI(fluidPage(
    
  # Application title
  titlePanel("Tabsets"),
  
  # Sidebar with controls to select the random distribution type
  # and number of observations to generate. Note the use of the
  # br() element to introduce extra vertical spacing
  sidebarLayout(
    sidebarPanel(
      radioButtons("dataset", "Data Set:",
                   c("Ovarian TCGA" = "ov",
                     "Pancreatic TCGA" = "panc",
                     "Prostate TCGA" = "pros",
                     "Liver TCGA" = "liv")),
      br(),
      
      textInput("gene", 
                  "Enter Gene", 
                 "121_at"),
	submitButton("Run")
    ),
    
    # Show a tabset that includes a plot, summary, and table view
    # of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs", 
        tabPanel("Volcano Plot", plotOutput("plot")), 
        tabPanel("Heatmap", plotOutput("summary")), 
        tabPanel("Table", tableOutput("table"))
      )
    )
  )
))
