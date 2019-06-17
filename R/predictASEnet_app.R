#' 
#' ASEnet_app
#' 
#' Shiny interface to run Predict.ASEnet function
#' 
#' NEED TO ADD ERROR HANDLING
#' 
#' Dependencies:
#' 
#'             ASEnet.R must be in same directory
#'             
#'             library(shiny)
#'             library(glmnet)
#'             
#'             

ui <- fluidPage(
  
  titlePanel("Predict.ASEnet"),
  
  p("Welcome to ASEnet, a machine learning R application designed to predict gene expression from nearby variants on the genome"),
  
  div(
    
    h4("Gene expression prediction"),
    
    p("This function uses a model created with the Train.ASEnet function to make predictions of the expression counts of either
      both alleles if the model is built at the chromosome-level, or of summed expression if the model is built at the genotype-level. 
      As input, a new chromosomal sequence of haplotypes must be provided")),
  
  selectInput("ASEmode", "Choose a prediction method: (must provide corresponding model and haplotypes)",
              list("Chromosome-level", 
                   "Genotype-level")
  ),
  
  fileInput("Model", "Choose model file",
            accept = ".rda"
  ),
  
  fileInput("newHaps", "Choose a file containing the new haplotypes",
            accept = ".rda"
  ),

  actionButton("do", "Predict"),
  textOutput("prediction"),
  p("\n"),p("\n")

)

server <- function(input, output) {
  
  options(shiny.maxRequestSize=6000*1024^2)
  
  observeEvent(input$do, {
    
    if (input$ASEmode == "Chromosome-level"){ASEmode <- 1}else{ASEmode <- 0}
    
    library(glmnet)
    source("./ASEnet.R")
    
    load(input$Model$datapath)
    load(input$newHaps$datapath)
    
    if(ASEmode == 1){Predict.ASEnet(ASEModel,newHapsASE, ASEmode = 1)}else{Predict.ASEnet(EModel,newHapsE, ASEmode = 0)}
    
  })

}

shinyApp(ui = ui, server = server)
