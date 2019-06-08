#' 
#' ASEnet
#' 
#' Machine learning algorithms to predict chromosome-level gene expression from allele specific rSNPs
#' 
#' Train.ASEnet: This function allows for the creation of chromosome-level gene expression 
#' predictive models based on the implementation of the lasso and elastic net penalised regression methodologies. 
#' The function uses chromosome-level rSNP data as factors in the prediction of coding SNP ASE information.
#' 
#' Arguments:
#' 
#'           gen_input - Rdata file, output of the Gen.input function present in the lettols GitHub repository.
#' 
#'           TSSwin - This represents the distance from each gene's TSS over which nearby variants will be selected,
#'                    either side of the transcript start site. Defaults to 500kb
#' 
#'           alpha - elastic net mixing parameter: 
#'                   a value of 1 adds the lasso penalty, while a value of 0 adds the ridge penalty, default of 0.5 for 
#'                   elastic net
#' 
#' 
#' Predict.ASEnet: Uses output model generated in Train.ASEnet, as well as chromosome-level rSNP data to make 
#' coding SNP ASE predictions
#' 
#' Arguments:
#' 
#'           ASEModel - Rdata file, output of the Train.Asenet function, containing the ASE model built
#'                          
#'           newHaps - Rdata file containing the haplotypes of both reference and alternative alleles of a new study 
#'                     individual. rSNPs must match those used for training (function recognises them by name)
#'                     
#' 
#' Dependencies:
#' 
#'             library(dplyr)
#'           
#'             library(glmnet)
#'             
#' 


Train.ASEnet <- function(gen_input,TSSwin = 5e+05, alpha = 0.5){
  
  
  # filter out infromation needed from Gene.Input output
  
  cat("\n Welcome to ASEnet, let's train a model for each gene using nearby variants and SNP ASE ...\n")
  
  rSNPs <- data.frame(ID=gen_input$leg$ID, end=gen_input$leg$end, gen_input$haps)
  
  ASE <- select(gen_input$ASE, ID, end, TSS, Gene, Ind, refCount, altCount)
  
  ASE <- ASE[order(ASE$TSS),]
  
  totReads <- gen_input$counts
  
  totReads <- totReads[order(totReads$Ind),]
  
  #' iterate through every unique gene available, taking TSSwin to look for nearby rSNPs, then train
  #' each ASE
  
  ASEModel <- list()
  
  i <- 1
  
  while (i <= length(unique(ASE$Gene))) {
    
    currGene <- as.character(unique(ASE$Gene)[[i]]) # current ASE gene
    
    # here the code is using the lowest of the TSS for each gene, might have to change this? (ask James)
    
    currTSS <- min(ASE$TSS[which(ASE$Gene == unique(ASE$Gene)[[i]])])
    
    cat("\n",round(i/length(unique(ASE$Gene))*100,2)," % completed >>")
    cat("\nCurrent ASE gene: ", currGene ,"with transcription start site: ",currTSS)
    
    tempRange <- seq(currTSS-TSSwin,currTSS+TSSwin)
    
    snpASE <- ASE[which(ASE$Gene == currGene),] # ASE data on that gene
    
    # nearby variants of the reference alleles on that gene
    nearVarsR <- t(rSNPs[which(rSNPs$end %in% tempRange), which(colnames(rSNPs) %in% snpASE$Ind),drop = FALSE]) 
    
    # nearby variants of the alternative alleles on that gene
    nearVarsA <- t(rSNPs[which(rSNPs$end %in% tempRange), which(colnames(rSNPs) %in% snpASE$Ind)+1,drop = FALSE])
    
    
    j <- 1
    
    
    while (j <= length(unique(snpASE$ID))) { #iterate through all ASE sites of a gene
      
      currASE <- snpASE[which(snpASE$ID == unique(snpASE$ID)[j]),] # ASE data on current ASE site
      
      currTot <- totReads[which(totReads$Ind %in% currASE$Ind),]
      
      
      # put total count data in correct format to normalise in next step
      x <-2
      while (x <= nrow(currTot)) {
        
        currTot[x-1,3] <- currTot[x,2] 
        
        x <- x +2
        
      }
      
      currTot <-currTot[which(currTot$All != 1),]
      
      colnames(currTot) <- c("Ind","refCount","altCount")
      
      #normalise
      currASE <- transform(currASE, totRef = currTot$refCount,totAlt= currTot$altCount,
                           refCountN = currASE$refCount / currTot$refCount, altCountN = currASE$altCount / currTot$altCount)
      
      cat("\nCurrent ASE site: ",unique(snpASE$ID)[j])
      
      
      # nearby variants of the reference alleles on that ASE site
      currVarsR <- as.matrix(nearVarsR[which(rownames(nearVarsR) %in% currASE$Ind),,drop = FALSE]) 
      
      # nearby variants of the alternative alleles on that ASE site
      currVarsA <- as.matrix(nearVarsA[which(rownames(nearVarsR) %in% currASE$Ind),,drop = FALSE])
      
      #join currVars and modify currASE for a unified model
      currVars <- rbind(currVarsR,currVarsA)
      currASET <- rbind(as.matrix(currASE$refCountN),as.matrix(currASE$altCountN))
      
      if (length(unique(currASET)) > 1){
        
        cat(" - ",length(currASET)," ASE data points used for model ")
        
        ASEModel[[currGene]][[unique(snpASE$ID)[j]]] <- glmnet(currVars,currASET, alpha=alpha)
        
      }
      
      else{
        
        cat("\n[ Warning - not enough ASE data points to create a model of this ASE site ]")
        
      }
      
      j <- j+1
      
    }
    
    i <- i+1
    
  }
  
  cat("\n\nASE model trained for this chromosome!\n\nSaving model file...")
  
  save(ASEModel, file = "./ASEModel.rda")
  
}

Predict.ASEnet <- function(ASEmodel,newHaps){
  
  # loop over models, looking for corresponding rSNPs to predict ASE
  
  allele0 <- list()
  allele1 <- list()
  
  i <- 1
  
  while (i<length(ASEModel)) {
    
    cat("\n",round(i/length(ASEModel)*100)," % completed >> Predicting ASE of gene ", names(ASEModel)[i])
    
    j <-1
    
    while (j<=length(ASEModel[[i]])) {
      
      cat("\nPredicting ASE counts of site:", names(ASEModel[[i]])[j])
      
      newVars <- t(newHaps[which(rownames(newHaps) %in% ASEModel[[i]][[j]][["beta"]]@Dimnames[[1]]),c(1,2)])
      
      # results for all different lambdas saved
      allele0[[names(ASEModel)[i]]][[names(ASEModel[[i]])[j]]] <- predict(ASEModel[[i]][[j]],t(newVars[1,]))
      allele1[[names(ASEModel)[i]]][[names(ASEModel[[i]])[j]]] <- predict(ASEModel[[i]][[j]],t(newVars[2,]))
      
      # results from last lambda used shown
      cat(" >> Prediction of allele 1:", allele0[[names(ASEModel)[i]]][[names(ASEModel[[i]])[j]]]
          [length(allele0[[names(ASEModel)[i]]][[names(ASEModel[[i]])[j]]])])
      cat(" >> Prediction of allele 2:", allele1[[names(ASEModel)[i]]][[names(ASEModel[[i]])[j]]]
          [length(allele1[[names(ASEModel)[i]]][[names(ASEModel[[i]])[j]]])])
      
      j <- j+1
      
    }
    
    i <- i+1
    
  }
  
  ASEpredict <-list(allele0 = allele0,allele1 = allele1)
  
  save(ASEpredict, file = "./predictASE.rda")
  
}