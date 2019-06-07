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
#'           ref/altModel - Rdata file, output of the Train.Asenet function, containing the reference alnd alternative 
#'                          ASE models built
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

  refModel <- list()
  altModel <- list()
    
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
    nearVarsR <- t(rSNPs[which(rSNPs$end %in% tempRange), which(colnames(rSNPs) %in% snpASE$Ind)]) 
    
    # nearby variants of the alternative alleles on that gene
    nearVarsA <- t(rSNPs[which(rSNPs$end %in% tempRange), which(colnames(rSNPs) %in% snpASE$Ind)+1]) 
    
    if (nrow(nearVarsR) > 1){
    
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
    currASE <- transform(currASE, totRef = currTot$refCount,totAlt= currTot$altCount, refCountN = currASE$refCount / currTot$refCount, altCountN = currASE$altCount / currTot$altCount)
    
    cat("\nCurrent ASE site: ",unique(snpASE$ID)[j])
    
    # nearby variants of the reference alleles on that ASE site
    currVarsR <- as.matrix(nearVarsR[which(rownames(nearVarsR) %in% currASE$Ind),]) 
    
    # nearby variants of the alternative alleles on that ASE site
    currVarsA <- as.matrix(nearVarsA[which(rownames(nearVarsR) %in% currASE$Ind),])
    
    if (length(unique(currASE[,6])) > 1){
      
      cat(" - ",length(currASE[,6])," ASE data points used for reference model ")
    
      refModel[[currGene]][[unique(snpASE$ID)[j]]] <- glmnet(currVarsR,currASEN[,10], alpha=alpha)
    
    }
    
    else{
      
      cat("\n[ Warning - not enough ASE data points to create a model for the reference allele of this ASE site ]")
      
    }
    
    if (length(unique(currASE[,7])) > 1){
      
      cat(" - ",length(currASE[,7])," ASE data points used for alternative model ")
      
      altModel[[currGene]][[unique(snpASE$ID)[j]]] <- glmnet(currVarsA,currASEN[,11], alpha=alpha)
      
    }
    
    else{
      
      cat("\n[ Warning - not enough ASE data points to create a model for the alternative allele of this ASE site ]")
      
    }
    
    j <- j+1
    
    }
    
    }
    
    else{
      
      cat(" [ Warning - not enough ASE data points to create a model for this gene ]")
      
    }
    
    i <- i+1
    
  }
  
  cat("\n\nReference and alternative models trained for this chromosome!\n\nSaving model files...")
  
  save(refModel, file = "./refModel.rda")
  save(altModel, file = "./altModel.rda")
  
}

Predict.ASEnet <- function(refModel,altModel,newHaps){
  
  library(glmnet)
  
  # loop over models, looking for corresponding rSNPs to predict ASE
  
  predictASERef <- list()
  predictASEAlt <- list()
  
  i <- 1
  
  while (i<length(refModel)) {
    
    cat("\n",round(i/length(refModel)*100)," % completed >> Predicting ASE of gene ", names(refModel)[i])
    
    j <-1
    
    while (j<=length(refModel[[i]])) {
      
      cat("\nPredicting reference ASE counts of site:", names(refModel[[i]])[j])
      
      newVars <- t(newHaps[which(rownames(newHaps) %in% refModel[[i]][[j]][["beta"]]@Dimnames[[1]]),1])
      
      predictASERef[[names(refModel)[i]]][[names(refModel[[i]])[j]]] <- predict(refModel[[i]][[j]],newVars, s=80)
      
      cat(" >> Prediction:", predictASERef[[names(refModel)[i]]][[names(refModel[[i]])[j]]])
      
      j <- j+1
      
    }
    
    i <- i+1
    
  }
  
  i <- 1
  
  while (i<length(altModel)) {
    
    cat("\n",round(i/length(altModel)*100)," % completed >> Predicting ASE of gene ", names(altModel)[i])
    
    j <-1
    
    while (j<=length(altModel[[i]])) {
      
      cat("\nPredicting alternative ASE counts of site:", names(altModel[[i]])[j])
      
      newVars <- t(newHaps[which(rownames(newHaps) %in% altModel[[i]][[j]][["beta"]]@Dimnames[[1]]),2])
      
        predictASEAlt[[names(altModel)[i]]][[names(altModel[[i]])[j]]] <- predict(altModel[[i]][[j]],newVars, s=80)
        
        cat(" >> Prediction:", predictASEAlt[[names(altModel)[i]]][[names(altModel[[i]])[j]]])
      
      j <- j+1
      
    }
    
    i <- i+1
    
  }

  save(predictASERef, file = "./predictASERef.rda")
  save(predictASEAlt, file = "./predictASEAlt.rda")
  
  
}
