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
#'             
#'             

ui <- fluidPage(
  
  titlePanel("Plot.ASEnet"),
  
  p("Welcome to ASEnet, a machine learning R application designed to predict gene expression from nearby variants on the genome"),
  
  div(
    
    h4("Accuracy plotting "),
    
    p("This function takes two cross validated ASEnet models and compares their accuracies across common trained expression locations.
      It outputs a plot with either the difference in mean cross validated errors between the models at each location or the actual
      error values per model, per location (measured in proportion of total expression in each site)")),
  
  fileInput("ASEModel", "Choose chromosome-level model file",
            accept = ".rda"
  ),
  
  fileInput("EModel", "Choose genotype-level model file",
            accept = ".rda"
  ),
  
  # fileInput("aseDat", "Choose aseDat file", (maybe not necessary to know the positions?)
  #           accept = ".rda"
  # ),
  
  selectInput("type", "Which plot do you want to be created?",
              list("MCVE Difference",
                   "MCVE Difference Ranking",
                   "MCVE")
  ),
  numericInput("thold", "Error threshold for MCVE Difference Ranking and MCVE plots: ", 0.001),

  actionButton("do", "Plot"),
  p("\n"),p("\n")

)

server <- function(input, output) {
  
  options(shiny.maxRequestSize=6000*1024^2)
  
  observeEvent(input$do, {
    
    if (input$type == "MCVE Difference"){type <- 1}
    else if(input$type == "MCVE Difference Ranking"){type <- 2}
    else if(input$type == "MCVE"){type <- 3}
    
    load(input$ASEModel$datapath)
    load(input$EModel$datapath)
    
    source("./ASEnet.R")
    
    Plot.ASEnet(ASEModel,EModel,type = type, thold = input$thold)
    
  })

}

shinyApp(ui = ui, server = server)
