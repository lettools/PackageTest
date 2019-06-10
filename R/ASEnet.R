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
#'           ASEmode - 1 for modelling of allele specific expression, 0 for genotype level modelling (like prediXcan)
#'           
#'           min - minimum number of ASE data points to create each SNP expression model 
#'                 (set to 10 to be able to perform 10 fold cv)
#'                 
#'           nfolds - number of cross validation folds (must match or exceed min value)
#' 
#' 
#' Predict.ASEnet: Uses output model generated in Train.ASEnet, as well as chromosome-level rSNP data to make 
#' coding SNP ASE predictions
#' 
#' Arguments:
#' 
#'           Model - Rdata file, output of the Train.Asenet function, containing the model built
#'                          
#'           newHaps - Rdata file containing the haplotypes of both reference and alternative alleles 
#'           (or combined genotype-level information) of a new study individual. rSNPs must match those used for training
#'           (function recognises them by name)
#'                     
#'           ASEmode - 1 for predicting of allele specific expression, 0 for genotype level modelling (like prediXcan)           
#'                     
#' 
#' Dependencies:
#' 
#'             library(dplyr)
#'           
#'             library(glmnet)
#'             
#' 


Train.ASEnet <- function(gen_input,TSSwin = 5e+05, alpha = 0.5, ASEmode = 1, min=10, nfolds = 10){
  
  
  # filter out infromation needed from Gene.Input output
  
  cat("\n Welcome to ASEnet, let's train a model for each gene using nearby variants and SNP corresponding expression values ...\n")
  
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
      
      # ASE mode
      if (ASEmode == 1){
        
        #join currVars and modify currASE for a unified model
        currVars <- rbind(currVarsR,currVarsA)
        currASET <- rbind(as.matrix(currASE$refCountN),as.matrix(currASE$altCountN))
        
        
      # PREDIXCAN mode
      }else{
        
        # sum both haplotypes for each rSNP and both expression levels in each ASE SNP
        currVars <-  currVarsR + currVarsA
        currASET <- as.matrix(currASE$refCountN) + as.matrix(currASE$altCountN)
        
      }
      
      if (length(unique(currASET)) >= min){
        
        cat(" - ",length(currASET)," data points used for model ")
        
        ASEModel[[currGene]][[unique(snpASE$ID)[j]]] <- cv.glmnet(currVars,currASET, alpha=alpha, nfolds=nfolds)
        
      }
      
      else{
        
        cat("\n[ Warning - not enough data points to create a model of this expression for this site ]")
        
      }
      
      j <- j+1
      
    }
    
    i <- i+1
    
  }
  
  if (ASEmode == 1){
  
    cat("\n\nASE model trained for this chromosome!\n\nSaving model file...")
  
    save(ASEModel, file = "./ASEModel.rda")
  
  }else{
    
    cat("\n\nExpression model trained for this chromosome!\n\nSaving model file...")
    
    EModel <- ASEModel
    
    save(EModel, file = "./EModel.rda")
    
    
  }
  
}

Predict.ASEnet <- function(Model,newHaps,ASEmode = 1){
  
  if (ASEmode == 1){
  
  # loop over models, looking for corresponding rSNPs to predict ASE
  
  allele0 <- list()
  allele1 <- list()
  
  i <- 1
  
  while (i<=length(Model)) {
    
    cat("\n",round(i/length(Model)*100)," % completed >> Predicting ASE of gene ", names(Model)[i])
    
    j <-1
    
    while (j<=length(Model[[i]])) {
      
      cat("\nPredicting ASE counts of site:", names(Model[[i]])[j])
      
      newVars <- t(newHaps[which(rownames(newHaps) %in% Model[[i]][[j]][["glmnet.fit"]][["beta"]]@Dimnames[[1]]),c(1,2)])
      
      # result for lamda.min saved
      allele0[[names(Model)[i]]][[names(Model[[i]])[j]]] <- predict(Model[[i]][[j]],t(newVars[1,]),s="lambda.min")
      allele1[[names(Model)[i]]][[names(Model[[i]])[j]]] <- predict(Model[[i]][[j]],t(newVars[2,]),s="lambda.min")
      
      cat(" >> Prediction of allele 1:", allele0[[names(Model)[i]]][[names(Model[[i]])[j]]])
      cat(" >> Prediction of allele 2:", allele1[[names(Model)[i]]][[names(Model[[i]])[j]]])
      
      j <- j+1
      
    }
    
    i <- i+1
    
  }
  
  ASEpredict <-list(allele0 = allele0,allele1 = allele1)
  
  save(ASEpredict, file = "./ASEpredict.rda")
  
  cat("\nASE predictions made! Please,check your source folder")
  
  }else{
    
    # loop over models, looking for corresponding rSNPs to predict genotype level expression
    
    Epredict <- list()
    
    i <- 1
    
    while (i<=length(Model)) {
      
      cat("\n",round(i/length(Model)*100)," % completed >> Predicting expression of gene ", names(Model)[i])
      
      j <-1
      
      while (j<=length(Model[[i]])) {
        
        cat("\nPredicting expression counts of site:", names(Model[[i]])[j])
        
        newVars <- t(newHaps[which(rownames(newHaps) %in% Model[[i]][[j]][["glmnet.fit"]][["beta"]]@Dimnames[[1]])])
        
        # result for lamda.min saved
        Epredict[[names(Model)[i]]][[names(Model[[i]])[j]]] <- predict(Model[[i]][[j]],t(newVars[1,]),s="lambda.min")
        
        cat(" >> Prediction:", allele0[[names(Model)[i]]][[names(Model[[i]])[j]]])
        
        j <- j+1
        
      }
      
      i <- i+1
      
    }
    
    save(Epredict, file = "./Epredict.rda")
    
    cat("\nPredictions made! Please,check your source folder")
 
  }
}

Plot.ASEnet <-function(Model){
  
  mcve <- list()
  
  i <- 1
  
  while (i<=length(Model)) {
    
    j <-1
    
    while (j<=length(Model[[i]])) {
      
      mcve[["Sites"]][[length(mcve[["Sites"]])+1]] <- names(Model[[i]])[j]
      
      mcve[["mcve"]][[length(mcve[["mcve"]])+1]] <- ASEModel[[i]][[j]]$cvm[which(ASEModel[[i]][[j]]$lambda == ASEModel[[i]][[j]]$lambda.min)]
      
      j <- j+1
      
    }
    
    i <- i+1
    
  }
  
  #have to make it prettier and documentation (maybe add display options)
  plot(mcve$mcve,ylab="MCVE", main="Mean cross validated error per expression site")
  
}