#' 
#' ASEnet_app
#' 
#' Shiny interface to run Train.ASEnet function
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
  
  titlePanel("Train.ASEnet"),
  
  p("Welcome to ASEnet, a machine learning R application designed to predict gene expression from nearby variants on the genome"),
  
  div(
    
    h4("Model Training"),
    
    p("In it's default mode, ASEnet uses allele specific regulatory variants to build a, cross validated, allele-specific (chromosome-level)
    expression model, but it is also capable of performing the modelling based on the sum of regulatory variant haplotypes and 
    total gene expression from both alleles (genotype-level)")),
  
  fileInput("aseDat", "Choose aseDat file",
            accept = ".rda"
  ),
  
  selectInput("ASEmode", "Choose a modelling method:",
              list("Chromosome-level", 
                   "Genotype-level")
  ),
  
  numericInput("TSSwin", "Select transcription start site window to look for regulatory variants (kBases): ", 500),
  numericInput("alpha", "Elastic net mixing parameter, from 0 (lasso) to 1 (ridge)", 0.5, min = 0, max = 1),
  numericInput("min", "Minimum number of expression data points to build the model at each position: ", 10, min=2),
  numericInput("nfolds", "Number of folds for cross-validation: ", 10),
  
  selectInput("halve", "Do you want to use only one allele, so accuracy is comparable with a genotype-level model? 
              (only applicable in chromosome-level modelling, do not use if modelling at a genotype level)",
              list("Yes",
                   "No")
  ),
  
  
  actionButton("do", "Train Model"),
  textOutput("modelTraining"),
  p("\n"),p("\n")

)

server <- function(input, output) {
  
  options(shiny.maxRequestSize=6000*1024^2)
  
  observeEvent(input$do, {
    
    if (input$ASEmode == "Chromosome-level"){ASEmode <- 1}else{ASEmode <- 0}
    
    if (input$halve == "Yes"){halve <- 1}else{halve <- 0}
    
    load(input$aseDat$datapath)
    
    source("./ASEnet.R")
    
    Train.ASEnet(aseDat, ASEmode = ASEmode, TSSwin = input$TSSwin*1000, alpha = input$alpha, min = input$min, nfolds = input$nfolds, halve = halve)
    
  })

}

shinyApp(ui = ui, server = server)
