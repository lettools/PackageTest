#' 
#' ASEnet_app
#' 
#' Shiny Web interface to run ASEnet functions
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
  
  titlePanel("ASEnet"),
  
  p("Welcome to ASEnet, an R application designed to predict gene expression from nearby variants on the genome"),
  
  div(
    
    h4("Model Training"),
    
    p("In it's default mode, ASEnet uses allele specific regulatory variants and expression to build the expression model,
      but it is also capable of performing modelling based on the sum of regulatory variant haplotypes and 
      total gene expression from both alleles")),
  
  fileInput("aseDat", "Choose aseDat file",
            accept = ".rda"
  ),
  
  selectInput("ASEmode", "Choose a modelling method:",
              list("Chromosome-level information", 
                   "Genotype-level information")
  ),
  
  numericInput("TSSwin", "Select transcription start site window (kBases): ", 500),
  numericInput("alpha", "Elastic net mixing parameter, from 0 (lasso) to 1 (ridge)", 0.5, min = 0, max = 1),
  numericInput("min", "Minimum number of expression data points to build the model at each position: ", 10, min=2),
  numericInput("nfolds", "Number of folds for cross-validation: ", 10),
  
  selectInput("halve", "Do you want to use only one allele, so accuracy is comparable with a genotype-level model? 
              (only applicable in chromosome level modelling)",
              list("Yes",
                   "No")
  ),
  
  
  actionButton("do", "Train Model"),
  textOutput("modelTraining")

)

server <- function(input, output) {
  
  options(shiny.maxRequestSize=6000*1024^2)
  
  observeEvent(input$do, {
    
    aseDat <- input$aseDat
    
    ASEmode <- input$ASEmode
    
    if (ASEmode == "Chromosome-level information"){ASEmode <- 1}else{ASEmode <- 0}
    
    TSSwin <- input$TSSwin
    
    alpha <- input$alpha
    
    min <- input$min
    
    nfolds <- input$nfolds
    
    halve <- input$halve
    
    if (halve == "Yes"){halve <- 1}else{halve <- 0}
    
    load(aseDat$datapath)
    
    source("./ASEnet.R")
    
    Train.ASEnet(aseDat, ASEmode = ASEmode, TSSwin = TSSwin, alpha = alpha, min = min, nfolds = nfolds, halve = halve)
    
  })

}

shinyApp(ui = ui, server = server)
