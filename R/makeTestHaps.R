

#' makeTestHaps
#'
#' Function to create test haplotypes from aseDat file. Creates one file of allele-specific haplotypes as found in aseDat
#' and another file containing combined haplotype data for both alleles of each individual
#'
#' Arguments
#'
#'          aseDat - aseDat file containing the haplotype information
#'

makeTestHaps <- function(aseDat, halve = 0) {
  testHapsASE <- aseDat$haps
  
  #extract and join haplotypes for genotype level prediction
  testHapsGen <- list()
  colnames <- list()
  
  i <- 1
  j <- 1
  
  while (i < length(testHapsASE)) {
    testHapsGen[[j]] <- testHapsASE[[i]] + testHapsASE[[i + 1]]
    colnames[[j]] <- colnames(testHapsASE)[[i]]
    i <- i + 2
    j <- j + 1
  }
  
  testHapsGen <- as.data.frame(testHapsGen)
  
  rownames(testHapsGen) <- rownames(testHapsASE)
  colnames(testHapsGen) <- colnames
  
  # rename and reorder haplotypes for prediction
  colnames(testHapsASE)[seq(0, length(testHapsASE), by = 2)] <-
    paste(colnames(testHapsASE)[seq(0, length(testHapsASE), by = 2)], " 1", sep = "")
  
  tempTHASE1 <-
    testHapsASE[which(grepl("\\s", colnames(testHapsASE))) - 1]
  tempTHASE2 <-
    testHapsASE[which(grepl("\\s", colnames(testHapsASE)))]
  
  if (halve == 0) {
    testHapsASE <- cbind(tempTHASE1, tempTHASE2)
    save(testHapsASE, file = "./testHapsASE.rda")
    
  } else if (halve == 1) {
    testHapsASE_halved <- tempTHASE1
    save(testHapsASE_halved, file = "./testHapsASE_halved.rda")
  }
  
  save(testHapsGen, file = "./testHapsGen.rda")
  
}