#' 
#' ASEnet_app
#' 
#' Shiny Web interface to run ASEnet functions
#' 
#' Dependencies:
#' 
#'             ASEnet.R must be in same directory
#'             
#'             library(shiny)
#'             
#'             

ui <- fluidPage(
  
  titlePanel("ASEnet"),
  
  sidebarLayout(
    sidebarPanel(
      p("Welcome to ASEnet, an R application designed to predict gene expression from nearby variants on the genome"),
                 
      p("In it's default mode, ASEnet uses allele specific regulatory variants and expression to build the expression model,
      but it is also capable of performing modelling based on the sum of regulatory variant haplotypes and 
      total gene expression from both alleles"),
    
      p("Data required for the modelling: output from gen.Input function from 
      the packagetest package from the lettols gitHub repository")),
    
    mainPanel(
      
      h4("Model Training"),
      
      fileInput("aseDat", "Choose aseDat file",
                accept = ".rda"
      ),
      
      
      selectInput("ASEmode", "Choose a modelling method:",
                  list("Chromosome-level information (allele-specific)", 
                       "Genotype-level information (non-allele-specific)")
      ),textOutput("mode"),
      
      
      actionButton("do", "Click Me"),
      textOutput("modelTraining")
      
      
    )
  )
)

server <- function(input, output) {
  
  options(shiny.maxRequestSize=6000*1024^2)
  
  observeEvent(input$do, {
  
    source("./ASEnet.R")
  
    aseDat <- input$aseDat
  
    mode <-output$mode
  
    if(mode == "Chromosome-level information (allele-specific)"){mode <- 1} else {mode <- 0}
  
    Train.ASEnet(aseDat,ASEmode = mode)
    
  })
    
}

shinyApp(ui = ui, server = server)
