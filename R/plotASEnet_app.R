#' 
#' ASEnet_app
#' 
#' Shiny interface to run Plot.ASEnet function
#' 
#' NEED TO ADD ERROR HANDLING
#' 
#' Dependencies:
#' 
#'             ASEnet.R must be in same directory
#'             
#'             library(shiny)
#'             library(dplyr)
#'             library(glmnet)
#'             
#'             

ui <- fluidPage(
  
  titlePanel("Plot.ASEnet"),
  
  p("Welcome to ASEnet, a machine learning R application designed to predict gene expression from nearby variants on the genome"),
  
  div(
    
    h4("Accuracy plotting "),
    
    p("This function takes two cross validated ASEnet models and compares their accuracies across common trained expression locations.
      It outputs a plot with either the difference in mean cross validated errors between the models at each location or the actual
      error values per model, per location (measured in proportion of total expression in each site). As an additional input, this 
      function requires the corresponding aseDat file to correlate expression location IDs to positions in the genome.")),
  
  fileInput("ASEModel", "Choose chromosome-level model file",
            accept = ".rda"
  ),
  
  fileInput("EModel", "Choose genotype-level model file",
            accept = ".rda"
  ),
  
  fileInput("aseDat", "Choose aseDat file",
            accept = ".rda"
  ),
  
  selectInput("type", "Which plot do you want to be created?",
              list("MCVE difference",
                   "MCVE")
  ),

  actionButton("do", "Plot"),
  p("\n"),p("\n")

)

server <- function(input, output) {
  
  options(shiny.maxRequestSize=6000*1024^2)
  
  observeEvent(input$do, {
    
    if (input$type == "MCVE difference"){type <- 1}else{type <- 2}
    
    load(input$ASEModel$datapath)
    load(input$EModel$datapath)
    load(input$aseDat$datapath)
    
    source("./ASEnet.R")
    
    Plot.ASEnet(ASEModel,EModel,aseDat,type = type)
    
  })

}

shinyApp(ui = ui, server = server)
