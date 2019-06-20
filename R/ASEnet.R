#'
#' ASEnet
#'
#' Machine learning algorithms to predict chromosome-level gene expression from allele specific rSNPs
#'
#' Train.ASEnet: This function allows for the creation of chromosome-level gene expression
#' predictive models based on the implementation of the lasso and elastic net penalised regression methodologies.
#' The function uses chromosome-level rSNP data as factors in the prediction of coding SNP ASE information.
#' Additionally, this function stores the expression values used for modelling at each position for later
#' prediction accuracy analyses
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
#'           halve - only use reference allele for modelling, to have same number of data points as genotype level
#'                   modelling, option needed for precise plotting as different number of data points will shift axis
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
#' Plot.ASEnet: Uses cross-validated models from train.ASEnet and plots mean cross validated error across the positions
#'              of all trained mutations
#'
#'Arguments:
#'
#'           EModel / ASEModel - Rdata files, output of the Train.Asenet function, containing the models built for
#'                              genotype level and chromosome level prediction accordingly
#'
#'           type - type 1 displays difference in accuracy between models across common mutations, type 2 shows a ranking
#'                  of the expression sites with a higher difference and type 3 shows actual mcve errors for each model
#'                  across common expression locations
#'
#'           thold - for types 2 and 3, mcve threshold of shown expression locations (type 3 displays any expression locations
#'                  for which any of the two models suprpasses the threshold)
#'
#'
#' Dependencies:
#'
#'             library(dplyr)
#'             library(glmnet)
#'


Train.ASEnet <-
  function(gen_input,
           ASEmode = 1,
           TSSwin = 5e+05,
           alpha = 0.5,
           min = 10,
           nfolds = 10,
           halve = 0) {
    cat(
      paste(
        "\nParameters selected: \nWindow ->",
        TSSwin,
        " \nAlpha -> ",
        alpha,
        " \nASEmode -> ",
        ASEmode,
        " \nMinumum number of expression data points -> ",
        min,
        "\nNumber of folds for cross-validation -> ",
        nfolds,
        "\nHalve data?: ",
        halve,
        "\n"
      )
    )
    
    # filter out infromation needed from Gene.Input output
    
    rSNPs <-
      data.frame(ID = gen_input$leg$ID,
                 end = gen_input$leg$end,
                 gen_input$haps)
    
    ASE <-
      select(gen_input$ASE, ID, end, TSS, Gene, Ind, refCount, altCount)
    
    ASE <- ASE[order(ASE$TSS),]
    
    totReads <- gen_input$counts
    
    totReads <- totReads[order(totReads$Ind),]
    
    #' iterate through every unique gene available, taking TSSwin to look for nearby rSNPs, then train
    #' each ASE
    
    ASEModel <- list()
    ASEexp <- list()
    
    i <- 1
    
    while (i <= length(unique(ASE$Gene))) {
      currGene <- as.character(unique(ASE$Gene)[[i]]) # current ASE gene
      
      # here the code is using the lowest of the TSS for each gene, might have to change this? (ask James)
      
      currTSS <-
        min(ASE$TSS[which(ASE$Gene == unique(ASE$Gene)[[i]])])
      
      cat("\n", round(i / length(unique(ASE$Gene)) * 100, 2), " % completed >>")
      cat("\nCurrent ASE gene: ",
          currGene ,
          "with transcription start site: ",
          currTSS)
      
      tempRange <- seq(currTSS - TSSwin, currTSS + TSSwin)
      
      snpASE <-
        ASE[which(ASE$Gene == currGene),] # ASE data on that gene
      
      # nearby variants of the reference alleles on that gene
      nearVarsR <-
        t(rSNPs[which(rSNPs$end %in% tempRange), which(colnames(rSNPs) %in% snpASE$Ind), drop = FALSE])
      
      # nearby variants of the alternative alleles on that gene
      nearVarsA <-
        t(rSNPs[which(rSNPs$end %in% tempRange), which(colnames(rSNPs) %in% snpASE$Ind) +
                  1, drop = FALSE])
      
      
      j <- 1
      
      
      while (j <= length(unique(snpASE$ID))) {
        #iterate through all ASE sites of a gene
        
        currASE <-
          snpASE[which(snpASE$ID == unique(snpASE$ID)[j]),] # ASE data on current ASE site
        
        currTot <- totReads[which(totReads$Ind %in% currASE$Ind),]
        
        
        # put total count data in correct format to normalise in next step
        x <- 2
        while (x <= nrow(currTot)) {
          currTot[x - 1, 3] <- currTot[x, 2]
          
          x <- x + 2
          
        }
        
        currTot <- currTot[which(currTot$All != 1),]
        
        colnames(currTot) <- c("Ind", "refCount", "altCount")
        
        #normalise
        currASE <-
          transform(
            currASE,
            totRef = currTot$refCount,
            totAlt = currTot$altCount,
            refCountN = currASE$refCount / currTot$refCount,
            altCountN = currASE$altCount / currTot$altCount
          )
        
        cat("\nCurrent ASE site: ", unique(snpASE$ID)[j])
        
        
        # nearby variants of the reference alleles on that ASE site
        currVarsR <-
          as.matrix(nearVarsR[which(rownames(nearVarsR) %in% currASE$Ind), , drop = FALSE])
        
        # nearby variants of the alternative alleles on that ASE site
        currVarsA <-
          as.matrix(nearVarsA[which(rownames(nearVarsR) %in% currASE$Ind), , drop = FALSE])
        
        # ASE mode
        if (ASEmode == 1) {
          if (halve == 0) {
            #join currVars and modify currASE for a unified model
            currVars <- rbind(currVarsR, currVarsA)
            currASET <-
              rbind(as.matrix(currASE$refCountN),
                    as.matrix(currASE$altCountN))
            
          } else{
            #only use half the data, to compare accuracy with genotype level model
            #(using reference allele only, might have to change this?)
            currVars <- currVarsR
            currASET <- as.matrix(currASE$refCountN)
            
            
          }
          
          
          # PREDIXCAN mode
        } else{
          # sum both haplotypes for each rSNP and both expression levels in each ASE SNP
          currVars <-  currVarsR + currVarsA
          currASET <-
            as.matrix(currASE$refCountN) + as.matrix(currASE$altCountN)
          
        }
        
        if (length(unique(currASET)) >= min) {
          cat(" - ", length(currASET), " data points used for model ")
          
          ASEModel[[currGene]][[unique(snpASE$ID)[j]]] <-
            cv.glmnet(currVars,
                      currASET,
                      alpha = alpha,
                      nfolds = nfolds)
          
          if(ASEmode ==1) {
            if (halve == 0) {
            ASEexp[[currGene]][[unique(snpASE$ID)[j]]] <-
              cbind(c(currASE$Ind, paste(currASE$Ind," 
                                         1", sep = "")), currASET)
            }
          }
          else{
            ASEexp[[currGene]][[unique(snpASE$ID)[j]]] <-
              cbind(currASE$Ind, currASET)
          }
  
        }
        
        else{
          cat(
            "\n[ Warning - not enough data points to create a model of this expression for this site ]"
          )
          
        }
        
        j <- j + 1
        
      }
      
      i <- i + 1
      
    }
    
    if (ASEmode == 1) {
      if (halve == 0) {
        cat("\n\nASE model trained for this chromosome!\n\nSaving files...")
        
        save(ASEModel, file = "./ASEModel.rda")
        save(ASEexp, file = "./ASEexp.rda")
        
      } else{
        cat("\n\nASE model trained for this chromosome!\n\nSaving files...")
        
        ASEModel_Halved <- ASEModel
        
        save(ASEModel, file = "./ASEModel_Halved.rda")
        save(ASEexp, file = "./ASEexp.rda")
        
        
      }
      
    } else{
      cat("\n\nExpression model trained for this chromosome!\n\nSaving files...")
      
      EModel <- ASEModel
      Eexp <- ASEexp
      
      save(EModel, file = "./EModel.rda")
      save(Eexp, file = "./Eexp.rda")
      
    }
    
  }


Predict.ASEnet <- function(Model, newHaps) {
  # loop over models, looking for corresponding rSNPs to predict genotype level expression
  
  predict <- list()
  
  i <- 1
  
  while (i <= length(Model)) {
    cat(
      "\n",
      round(i / length(Model) * 100),
      " % completed >> Predicting expression of gene ",
      names(Model)[i]
    )
    
    j <- 1
    
    while (j <= length(Model[[i]])) {
      cat("\nPredicting expression counts of site:",
          names(Model[[i]])[j])
      
      newVars <-
        t(newHaps[which(rownames(newHaps) %in% Model[[i]][[j]][["glmnet.fit"]][["beta"]]@Dimnames[[1]]),])
      
      # result for lamda.min saved
      predict[[names(Model)[i]]][[names(Model[[i]])[j]]] <-
        cbind(rownames(newVars),
              predict(Model[[i]][[j]], newVars, s = "lambda.min"))
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  save(predict, file = "./predict.rda")
  
  cat("\n\nPredictions made! Please,check your source folder")
  
}

Plot.mcve.ASEnet <- function(ASEModel,
                        EModel,
                        aseDat,
                        type = 1,
                        thold = 0.001) {
  # set margins for better visibility
  par(mar = c(6.2, 6, 2, 2) + 0.1, mgp = c(4.5, 1, 0))
  
  # retrieve ASE mean cross validated error values
  ASEmcve <- list()
  
  i <- 1
  
  while (i <= length(ASEModel)) {
    j <- 1
    
    while (j <= length(ASEModel[[i]])) {
      ASEmcve[["Sites"]][[length(ASEmcve[["Sites"]]) + 1]] <-
        names(ASEModel[[i]])[j]
      
      ASEmcve[["mcve"]][[length(ASEmcve[["mcve"]]) + 1]] <-
        ASEModel[[i]][[j]]$cvm[which(ASEModel[[i]][[j]]$lambda == ASEModel[[i]][[j]]$lambda.min)]
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  # retrieve genotype level expression mean cross validated error values
  Emcve <- list()
  
  i <- 1
  
  while (i <= length(EModel)) {
    j <- 1
    
    while (j <= length(EModel[[i]])) {
      Emcve[["Sites"]][[length(Emcve[["Sites"]]) + 1]] <-
        names(EModel[[i]])[j]
      
      Emcve[["mcve"]][[length(Emcve[["mcve"]]) + 1]] <-
        EModel[[i]][[j]]$cvm[which(EModel[[i]][[j]]$lambda == EModel[[i]][[j]]$lambda.min)]
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  # retrie expresssion position from site name
  
  ASEmcve[["Locations"]] <-
    aseDat$ASE$end[which(ASEmcve$Sites %in% aseDat$ASE$ID)]
  ASEmcve$Locations <- ASEmcve$Locations[order(ASEmcve$Locations)]
  
  Emcve[["Locations"]] <-
    aseDat$ASE$end[which(Emcve$Sites %in% aseDat$ASE$ID)]
  Emcve$Locations <- Emcve$Locations[order(Emcve$Locations)]
  
  #get common expression sites for plot
  ASEmcveP <- list()
  ASEmcveP[["Sites"]] <-
    ASEmcve$Sites[which(ASEmcve$Sites %in% Emcve$Sites)]
  ASEmcveP[["mcve"]] <-
    ASEmcve$mcve[which(ASEmcve$Sites %in% Emcve$Sites)]
  ASEmcveP[["Locations"]] <-
    ASEmcve$Locations[which(ASEmcve$Sites %in% Emcve$Sites)]
  
  EmcveP <- list()
  EmcveP[["Sites"]] <-
    Emcve$Sites[which(Emcve$Sites %in% ASEmcve$Sites)]
  EmcveP[["mcve"]] <-
    Emcve$mcve[which(Emcve$Sites %in% ASEmcve$Sites)]
  EmcveP[["Locations"]] <-
    Emcve$Locations[which(Emcve$Sites %in% ASEmcve$Sites)]
  
  
  # actual plotting
  
  # mcve difference
  if (type == 1) {
    plot(
      EmcveP$mcve - ASEmcveP$mcve,
      type = "l",
      xlab = "Common expression sites",
      ylab = "Emodel mcve - ASEModel mcve"
    )
    
    cat("\n\nPlot created, please check the plot window in RStudio")
    
  # mcve ranked with a threshold
  } else if (type == 2) {
    errDif <- list()
    errDif <- EmcveP$mcve - ASEmcveP$mcve
    
    errDif <- list(difs = errDif, sites = ASEmcveP$Sites)
    
    errDif$sites <- errDif$sites[order(-errDif$difs)]
    errDif$difs <- errDif$difs[order(-errDif$difs)]
    
    errDif$sites <- errDif$sites[which(abs(errDif$difs) > thold)]
    errDif$difs <- errDif$difs[which(abs(errDif$difs) > thold)]
    
    
    barplot(
      errDif$difs,
      names.arg = errDif$sites,
      main = "Common expression sites",
      ylab = paste(
        "Emodel mcve - ASEModel mcve (ranked, threshold of :",
        thold,
        ")"
        , sep = ""),
      las = 2,
      cex.names = 1,
      xpd = FALSE,
      beside = TRUE,
      ylim = range(pretty(c(0, errDif$difs)))
    )
    
    cat("\n\nPlot created, please check the plot window in RStudio")
    
    
  # actual mcve values with a threshold
  } else if (type == 3) {
    blue <- rgb(0, 0, 1, alpha = 0.5)
    red <- rgb(1, 0, 0, alpha = 0.5)
    
    commThold <-
      unique(c(which(abs(ASEmcveP$mcve) > thold), which(abs(EmcveP$mcve) > thold)))
    commThold <- commThold[order(commThold)]
    
    ASEmcveP$Sites <- ASEmcveP$Sites[commThold]
    ASEmcveP$mcve <- ASEmcveP$mcve[commThold]
    EmcveP$Sites <- EmcveP$Sites[commThold]
    EmcveP$mcve <- EmcveP$mcve[commThold]
    
    maxmcve <- max(max(ASEmcveP$mcve), max(EmcveP$mcve))
    
    barplot(
      EmcveP$mcve,
      names.arg =  ASEmcveP$Sites,
      main = "Common expression sites",
      ylab = paste("mcve (threshold of :", thold, ")", sep = ""),
      las = 2,
      cex.names = 1,
      col = red,
      ylim = range(pretty(c(0, maxmcve)))
    )
    legend(
      "topleft",
      legend = c("Genotype-level", "Chromosome-level"),
      fill = c(red, blue),
      bty = "n",
      cex = 0.75
    )
    
    par(new = TRUE)
    barplot(
      ASEmcveP$mcve,
      col = blue,
      axes = FALSE,
      ylim = range(pretty(c(0, maxmcve)))
    )
    
    cat("\n\nPlot created, please check the plot window in RStudio")
    
  }
  
}

Plot.R2.ASEnet <- function(predict, exp) {
  i <- 1
  
  # retrieve predictions of known expressions to calculate R2
  while (i <= length(exp)) {
    
    j <- 1
    
    while (j <= length(exp[[i]])) {
      predict[[i]][[j]] <-
        predict[[i]][[j]][which(predict[[i]][[j]][, 1] %in% exp[[i]][[j]][, 1])]
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  # calaculate R squared
  R2 <-  list()
  
  while (i <= length(exp)) {
    
    j <- 1
    
    while (j <= length(exp[[i]])) {
      
      R2[[names(exp)[i]]][[names(exp[[i]])[j]]] <-
        cor(predict[[i]][[j]],exp[[i]][[j]]) ^ 2
  
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  # finish plotting
  
  cat("Plot created, please chack the plot window in RStudio")
  
}