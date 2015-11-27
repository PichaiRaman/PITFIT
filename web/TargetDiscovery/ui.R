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
                     "Pancreatic TCGA" = "paad",
                     "Prostate TCGA" = "prad")),
      br(),
      
      selectizeInput('gene', 'Enter Feature', choices = 'FOXM1_Exp', multiple=FALSE),

      br(),
      
      textInput("pval", 
                  "Adjusted P-Value Cutoff", 
                 ".25"),
      br(),

      textInput("logFC", 
                  "Log FC Cutoff", 
                 "1"),
	submitButton("Run"),
br(),
downloadButton('downloadData', 'Download')
	
    ),

    
    # Show a tabset that includes a plot, summary, and table view
    # of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs", 
        tabPanel("Table", dataTableOutput(outputId="table")), 
       tabPanel("Volcano Plot", plotOutput("plot")), 
        tabPanel("Heatmap", plotOutput("summary")) 
      )
    )
  )
))
