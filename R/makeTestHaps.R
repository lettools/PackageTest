
#' makeTestHaps
#' 
#' Function to create test haplotypes from aseDat file. Creates one file of allele-specific haplotypes as found in aseDat
#' and another file containing combined haplotype data for both alleles of each individual
#' 
#' Arguments
#' 
#'          aseDat - aseDat file containing the haplotype information
#'                    

makeTestHaps <- function(aseDat){
  
  testHapsASE <- aseDat$haps
  testHapsGen <- list()
  colnames <- list()
  
  i <-1
  j <-1
  
  while (i < length(testHapsASE)) {
    testHapsGen[[j]] <- testHapsASE[[i]]+testHapsASE[[i+1]]
    colnames[[j]] <- colnames(testHapsASE)[[i]]
    i <- i+2
    j <- j+1
  }
  
  testHapsGen <- as.data.frame(testHapsGen)
  
  rownames(testHapsGen) <- rownames(testHapsASE)
  colnames(testHapsGen) <- colnames
  
  save(testHapsASE, file = "./testHapsASE.rda")
  save(testHapsGen, file = "./testHapsGen.rda")
  
}