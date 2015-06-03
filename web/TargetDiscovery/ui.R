library(shiny)

# Define UI for random distribution application 
shinyUI(fluidPage(
    
  # Application title
  titlePanel("Target Discovery"),
  
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

      br(),
      
      textInput("pval", 
                  "Adjusted P-Value Cutoff", 
                 ".25"),
      textInput("logFC", 
                  "Log FC Cutoff", 
                 "1"),
	submitButton("Run")
	
    ),
    
    # Show a tabset that includes a plot, summary, and table view
    # of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs", 
        tabPanel("Table", tableOutput("table")), 
       tabPanel("Volcano Plot", plotOutput("plot")), 
        tabPanel("Heatmap", plotOutput("summary")) 
      )
    )
  )
))
