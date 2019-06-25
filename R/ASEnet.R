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
#'           aseDat - Rdata file, output of the Gen.input function present in the lettols GitHub repository.
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
#'                   modelling
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
#'           lambda - lambda reduction parameter to be used. mincverr selects lambda producing a minimum cross
#'                    validated error during model creation, se selects largest value of lambda such that error is within
#'                    1 standard error of the minimum, lowest selects lowest value of lambda used in model creation
#'
#'
#'Plot.mcve.ASEnet: Uses cross-validated models from train.ASEnet and plots mean cross validated error across the positions
#'                  of all trained mutations
#'
#'Arguments:
#'
#'           GenModel / ASEModel - Rdata files, output of the Train.Asenet function, containing the models built for
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
#'Plot.R2.ASEnet: Plots observed expression values against predicted values. Also calculates and plots the value of r squared
#'                per expression site
#'
#'Arguments:
#'
#'           model - ASE or genotype level model file, containing expression observed at each modelled site
#'           
#'           predict - predictions of gene expression values, outputted from the Predict.ASEnet function
#'
#'           type - type 1 plots the R2 value across gene expressionsites while type 2 shows an observed vs 
#'                  predicted scatter plot
#'                  
#'           sitename - (for type 2 plot) name of the expression site for wich to see observed vs predicted plot
#'
#'
#' Dependencies:
#'
#'             library(dplyr)
#'             library(glmnet)
#'


Train.ASEnet <-
  function(aseDat,
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
      data.frame(ID = aseDat$leg$ID,
                 end = aseDat$leg$end,
                 aseDat$haps)
    
    ASE <-
      select(aseDat$ASE, ID, end, TSS, Gene, Ind, refCount, altCount)
    
    ASE <- ASE[order(ASE$TSS), ]
    
    totReads <- aseDat$counts
    
    totReads <- totReads[order(totReads$Ind), ]
    
    #' iterate through every unique gene available, taking TSSwin to look for nearby rSNPs, then train
    #' each ASE
    
    ASEModel <- list()
    
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
        ASE[which(ASE$Gene == currGene), ] # ASE data on that gene
      
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
          snpASE[which(snpASE$ID == unique(snpASE$ID)[j]), ] # ASE data on current ASE site
        
        currTot <- totReads[which(totReads$Ind %in% currASE$Ind), ]
        
        
        # put total count data in correct format to normalise in next step
        
        x <- 2
        while (x <= nrow(currTot)) {
          currTot[x - 1, 3] <- currTot[x, 2]
          
          x <- x + 2
          
        }
        
        currTot <- currTot[which(currTot$All != 1), ]
        
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
        
        # glmnet implementation, locations and expression saving for plotting / R2 calculation
        
        if (length(unique(currASET)) >= min) {
          cat(" - ", length(currASET), " data points used for model ")
          
          ASEModel[[currGene]][[unique(snpASE$ID)[j]]][["model"]] <-
            cv.glmnet(currVars,
                      currASET,
                      alpha = alpha,
                      nfolds = nfolds)
          
          if (ASEmode == 1 && halve == 0) {
            ASEModel[[currGene]][[unique(snpASE$ID)[j]]][["expression"]] <-
              cbind(c(currASE$Ind, paste(currASE$Ind, "1", sep = " ")), currASET)
          }
          else{
            ASEModel[[currGene]][[unique(snpASE$ID)[j]]][["expression"]] <-
              cbind(currASE$Ind, currASET)
          }
          
            
          ASEModel[[currGene]][[unique(snpASE$ID)[j]]][["location"]] <- currASE$end[1]
          
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
        
      } else{
        cat("\n\nASE model trained for this chromosome!\n\nSaving files...")
        
        ASEModel_Halved <- ASEModel
        
        save(ASEModel_Halved, file = "./ASEModel_Halved.rda")
        
      }
      
    } else{
      cat(
        "\n\nGenotype level expression model trained for this chromosome!\n\nSaving files..."
      )
      
      GenModel <- ASEModel
      
      save(GenModel, file = "./GenModel.rda")
      
    }
    
  }


Predict.ASEnet <-
  function(Model,
           newHaps,
           ASEmode,
           lambda = "mincverr") {
    # loop over models, looking for corresponding rSNPs to predict genotype level expression
    
    if (is.numeric(ASEmode) == TRUE) {
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
            t(newHaps[which(rownames(newHaps) %in% Model[[i]][[j]][["model"]][["glmnet.fit"]][["beta"]]@Dimnames[[1]]), ])
          
          # selection of lambda
          
          if (lambda == "mincverr") {
            lambda <- "lambda.min"
          }
          else if (lambda == "se") {
            lambda <- "lambda.1se"
          }
          else if (lambda == "lowest") {
            lambda <-
              Model[[i]][[j]][["model"]]$lambda[length(Model[[i]][[j]][["model"]]$lambda)]
          }
          
          # predictions
          
          predict[[names(Model)[i]]][[names(Model[[i]])[j]]] <-
            cbind(rownames(newVars),
                  predict(Model[[i]][[j]][["model"]], newVars, s = lambda))
          
          j <- j + 1
          
        }
        
        i <- i + 1
        
      }
      
      if (ASEmode == 1) {
        ASEpredict <- predict
        
        save(ASEpredict, file = "./ASEpredict.rda")
        
        cat("\n\nASE predictions made! Please,check your source folder")
        
      } else{
        Genpredict <- predict
        
        save(Genpredict, file = "./Genpredict.rda")
        
        cat(
          "\n\nGenotype level expression predictions made! Please,check your source folder"
        )
        
      }
      
    }
    
  }

Plot.mcve.ASEnet <- function(GenModel,
                             ASEModel,
                             type = 1,
                             thold = 0.001) {
  # set margins for better visibility
  
  par(mar = c(6.2, 6, 2, 2) + 0.1, mgp = c(4.5, 1, 0))
  
  # retrieve genotype level expression mean cross validated error values
  
  Genmcve <- list()
  
  i <- 1
  
  while (i <= length(GenModel)) {
    j <- 1
    
    while (j <= length(GenModel[[i]])) {
      Genmcve[["sites"]][[length(Genmcve[["sites"]]) + 1]] <-
        names(GenModel[[i]])[j]
      
      Genmcve[["locations"]][[length(Genmcve[["sites"]]) + 1]] <-
        GenModel[[i]][[j]][["location"]]
      
      Genmcve[["mcve"]][[length(Genmcve[["mcve"]]) + 1]] <-
        GenModel[[i]][[j]][["model"]]$cvm[which(GenModel[[i]][[j]][["model"]]$lambda == GenModel[[i]][[j]][["model"]]$lambda.min)]
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  # retrieve ASE mean cross validated error values (these two steps could be put in a loop)
  
  ASEmcve <- list()
  
  i <- 1
  
  while (i <= length(ASEModel)) {
    j <- 1
    
    while (j <= length(ASEModel[[i]])) {
      ASEmcve[["sites"]][[length(ASEmcve[["sites"]]) + 1]] <-
        names(ASEModel[[i]])[j]
      
      ASEmcve[["locations"]][[length(ASEmcve[["sites"]]) + 1]] <-
        ASEModel[[i]][[j]][["location"]]
      
      ASEmcve[["mcve"]][[length(ASEmcve[["mcve"]]) + 1]] <-
        ASEModel[[i]][[j]][["model"]]$cvm[which(ASEModel[[i]][[j]][["model"]]$lambda == ASEModel[[i]][[j]][["model"]]$lambda.min)]
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  #get common expression sites for plot
  
  GenmcveC <- list()
  GenmcveC[["sites"]] <-
    Genmcve$sites[which(Genmcve$sites %in% ASEmcve$sites)]
  GenmcveC[["locations"]] <-
    Genmcve$locations[which(Genmcve$sites %in% ASEmcve$sites)]
  GenmcveC[["mcve"]] <-
    Genmcve$mcve[which(Genmcve$sites %in% ASEmcve$sites)]
  
  ASEmcveC <- list()
  ASEmcveC[["sites"]] <-
    ASEmcve$sites[which(ASEmcve$sites %in% Genmcve$sites)]
  ASEmcveC[["locations"]] <-
    ASEmcve$locations[which(ASEmcve$sites %in% Genmcve$sites)]
  ASEmcveC[["mcve"]] <-
    ASEmcve$mcve[which(ASEmcve$sites %in% Genmcve$sites)]
  
  
  
  # actual plotting
  
  # mcve difference
  
  if (type == 1) {
    plot(
      ASEmcveC$locations,
      GenmcveC$mcve - ASEmcveC$mcve,
      type = "p",
      xlab = "Chromosome location",
      ylab = "GenModel mcve - ASEModel mcve"
    )
    
    cat("\n\nPlot created, please check the plot window in RStudio")
    
    # mcve ranked with a threshold
    
  } else if (type == 2) {
    errDif <- list()
    errDif <- GenmcveC$mcve - ASEmcveC$mcve
    
    errDif <-
      list(
        difs = errDif,
        sites = ASEmcveC$sites,
        locations = ASEmcveC$locations
      )
    
    errDif$sites <- errDif$sites[order(-errDif$difs)]
    errDif$locations <- errDif$locations[order(-errDif$difs)]
    errDif$difs <- errDif$difs[order(-errDif$difs)]
    
    errDif$sites <- errDif$sites[which(abs(errDif$difs) > thold)]
    errDif$locations <-
      errDif$locations[which(abs(errDif$difs) > thold)]
    errDif$difs <- errDif$difs[which(abs(errDif$difs) > thold)]
    
    
    barplot(
      errDif$difs,
      names.arg = paste(errDif$sites, errDif$locations, sep = "\n"),
      main = "Common expression sites",
      ylab = paste(
        "GenModel mcve - ASEModel mcve (ranked, threshold of :",
        thold,
        ")"
        ,
        sep = ""
      ),
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
      unique(c(which(abs(ASEmcveC$mcve) > thold), which(abs(GenmcveC$mcve) > thold)))
    commThold <- commThold[order(commThold)]
    
    ASEmcveC$sites <- ASEmcveC$sites[commThold]
    ASEmcveC$locations <- ASEmcveC$locations[commThold]
    ASEmcveC$mcve <- ASEmcveC$mcve[commThold]
    
    GenmcveC$sites <- GenmcveC$sites[commThold]
    GenmcveC$locations <- GenmcveC$locations[commThold]
    GenmcveC$mcve <- GenmcveC$mcve[commThold]
    
    maxmcve <- max(max(ASEmcveC$mcve), max(GenmcveC$mcve))
    
    barplot(
      GenmcveC$mcve,
      names.arg =  paste(ASEmcveC$sites, ASEmcveC$locations, sep = "\n"),
      main = "Common expression sites",
      ylab = paste("mcve (threshold of :", thold, ")", sep = ""),
      las = 2,
      cex.names = 0.75,
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
      ASEmcveC$mcve,
      col = blue,
      axes = FALSE,
      ylim = range(pretty(c(0, maxmcve)))
    )
    
    cat("\n\nPlot created, please check the plot window in RStudio")
    
  }
  
}

Plot.R2.ASEnet <- function(Model, predict, type = 1, sitename = "rs2845405") {
  i <- 1
  
  # retrieve predictions of known expressions to calculate R2
  
  while (i <= length(Model)) {
    j <- 1
    
    while (j <= length(Model[[i]])) {
      predict[[i]][[j]] <-
        predict[[i]][[j]][which(predict[[i]][[j]][, 1] %in% Model[[i]][[j]][["expression"]][, 1]),]
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  # calaculate R squared

  R2 <-  list()
  observed <- list()
  predicted <- list()
  
  i <- 1
  
  while (i <= length(Model)) {
    j <- 1
    
    while (j <= length(Model[[i]])) {
      R2[["sites"]][[length(R2[["sites"]]) + 1]] <-
        names(Model[[i]])[j]
      
      R2[["locations"]][[length(R2[["locations"]]) + 1]] <-
        Model[[i]][[j]]$location
      
      R2[["people"]][[length(R2[["people"]]) + 1]] <-
        Model[[i]][[j]]$model$glmnet.fit$nobs
      
      R2[["values"]][[length(R2[["values"]]) + 1]] <-
        summary(lm(as.numeric(predict[[i]][[j]][, 2]) ~ as.numeric(Model[[i]][[j]][["expression"]][, 2])))$r.squared
      
      # retrieve observed abd predicted values
      
      observed[names(Model[[i]])[j]] <-
        list(as.numeric(Model[[i]][[j]][["expression"]][, 2]))
      
      predicted[names(Model[[i]])[j]] <-
        list(as.numeric(predict[[i]][[j]][, 2]))
      
      
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  if (type == 1) {
    
    par()
    dev.off()
    
    # set bins and colour pallete for gradient depending on number of individuals
    
    bins <- cut(R2[["people"]],10,include.lowest=TRUE)
    colfunc <- colorRampPalette(c("grey", "black"))
    pallete <- colfunc(10)
  
    plotcols <- pallete[bins]
    
    par(mar=c(5.1, 4.1, 4.1, 10), xpd=TRUE)
    
    plot(R2[["locations"]],R2[["values"]], xlab = "expression locations", ylab = "R2", col = plotcols, pch = 16, bty = "n")
    
    legend(
      "topright",
      inset=c(-0.36,0),
      legend = levels(bins),
      pch = 16,
      col = pallete,
      bty = "n",
      cex = 1
    )
    
    
  } else if (type == 2) {
    
    par()
    dev.off()
    
    plot(array(as.numeric(unlist(observed[sitename]))), array(as.numeric(unlist(predicted[sitename]))), xlab = "observed", 
         ylab = "predicted", main = sitename)
    
  }
  
  cat("Plot created, please check the plot window in RStudio")
  
}