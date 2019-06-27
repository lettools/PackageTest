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
#'                     (just to save file with matching name)
#'
#'           halve - is the ASE model halved (just to save file with matching name), only needed if ASEMode = 1
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
#'           test - "R2" for r squared testing and "SPM" for spearman testing
#'
#'           type - type 1 plots the R2 value across gene expressionsites, type 2 shows an observed vs
#'                  predicted scatter plot, type 3 plots R2 / SPM values of a model against another model
#'                  (secmodel and secpredict needed)
#'
#'           sModel/spredict - model and predict files for secondary model for tyoe 3 plotting
#'
#'
#'
#' Dependencies:
#'
#'             library(tidyverse)
#'             library(reshape2)
#'             library(glmnet)
#'
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
          
          
          ASEModel[[currGene]][[unique(snpASE$ID)[j]]][["location"]] <-
            currASE$end[1]
          
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
           halve,
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
        if (halve == 0) {
          ASEpredict <- predict
          
          save(ASEpredict, file = "./ASEpredict.rda")
          
          cat("\n\nASE predictions made! Please,check your source folder")
          
        } else{
          ASEpredict_Halved <- predict
          
          save(ASEpredict_Halved, file = "./ASEpredict.rda")
          
          cat("\n\nASE predictions made! Please,check your source folder")
          
        }
        
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
      
      Genmcve[["people"]][[length(Genmcve[["sites"]]) + 1]] <-
        GenModel[[i]][[j]][["model"]]$glmnet.fit$nobs
      
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
      
      ASEmcve[["people"]][[length(ASEmcve[["sites"]]) + 1]] <-
        ASEModel[[i]][[j]][["model"]]$glmnet.fit$nobs
      
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
  GenmcveC[["people"]] <-
    Genmcve$people[which(Genmcve$sites %in% ASEmcve$sites)]
  GenmcveC[["mcve"]] <-
    Genmcve$mcve[which(Genmcve$sites %in% ASEmcve$sites)]
  
  ASEmcveC <- list()
  ASEmcveC[["sites"]] <-
    ASEmcve$sites[which(ASEmcve$sites %in% Genmcve$sites)]
  ASEmcveC[["locations"]] <-
    ASEmcve$locations[which(ASEmcve$sites %in% Genmcve$sites)]
  ASEmcveC[["people"]] <-
    ASEmcve$people[which(ASEmcve$sites %in% Genmcve$sites)]
  ASEmcveC[["mcve"]] <-
    ASEmcve$mcve[which(ASEmcve$sites %in% Genmcve$sites)]
  
  # actual plotting
  
  # mcve difference
  
  if (type == 1) {
    data = data.frame(
      locations = ASEmcveC$locations,
      Gen_ASE_Model_MCVE_differences = GenmcveC$mcve - ASEmcveC$mcve,
      peopleGen = GenmcveC$people,
      peopleASE = ASEmcveC$people
    )
    
    p <- ggplot(data = data) +
      geom_jitter(aes(locations, Gen_ASE_Model_MCVE_differences))
    
    print(p)
    
    cat("\nPlot created, please check the plot window in RStudio\n ")
    
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
    
    errDif$sites <- factor(errDif$sites, levels = errDif$sites)
    
    errDif$sites <- errDif$sites[which(abs(errDif$difs) > thold)]
    errDif$locations <-
      errDif$locations[which(abs(errDif$difs) > thold)]
    errDif$difs <- errDif$difs[which(abs(errDif$difs) > thold)]
    
    p <-
      ggplot(data = data.frame(errDif), aes(sites, difs)) + geom_col()
    
    print(p)
    
    cat("\nPlot created, please check the plot window in RStudio\n ")
    
    
    # actual mcve values with a threshold
    
  } else if (type == 3) {
    commThold <-
      unique(c(which(abs(ASEmcveC$mcve) > thold), which(abs(GenmcveC$mcve) > thold)))
    commThold <- commThold[order(commThold)]
    
    ASEmcveC$sites <- ASEmcveC$sites[commThold]
    ASEmcveC$locations <- ASEmcveC$locations[commThold]
    ASEmcveC$mcve <- ASEmcveC$mcve[commThold]
    
    ASEmcveC$sites <-
      factor(ASEmcveC$sites, levels = ASEmcveC$sites)
    
    GenmcveC$sites <- GenmcveC$sites[commThold]
    GenmcveC$locations <- GenmcveC$locations[commThold]
    GenmcveC$mcve <- GenmcveC$mcve[commThold]
    
    maxmcve <- max(max(ASEmcveC$mcve), max(GenmcveC$mcve))
    
    data = data.frame(
      sites = ASEmcveC$sites,
      locations = ASEmcveC$locations,
      GenMCVE = GenmcveC$mcve,
      ASEMCVE = ASEmcveC$mcve
    )
    
    data <-
      melt(data[, c('sites' , 'GenMCVE', "ASEMCVE")], id.vars = 1)
    
    colnames(data)[2] <- "Models"
    colnames(data)[3] <- "MCVE"
    
    p <- ggplot(data
                , aes(sites, MCVE)) +
      geom_bar(aes(fill = Models), stat = "identity", position = "dodge") +
      theme(axis.text.x = element_text(angle = 90)) +
      labs(x = paste("Common sites with a MCVE lower threshold of ", thold, sep = ""))
    
    print(p)
    
    cat("\nPlot created, please check the plot window in RStudio\n ")
    
  }
  
}

Plot.R2.ASEnet <- function(Model,
                           predict,
                           type = 1,
                           test = "R2",
                           sModel,
                           spredict) {
  i <- 1
  
  # retrieve predictions of known expressions to calculate R2
  
  while (i <= length(Model)) {
    j <- 1
    
    while (j <= length(Model[[i]])) {
      predict[[i]][[j]] <-
        predict[[i]][[j]][which(predict[[i]][[j]][, 1] %in% Model[[i]][[j]][["expression"]][, 1]), ]
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  # calaculate R squared
  
  R2 <-  list()
  predicted <- list()
  observed <- list()
  
  
  i <- 1
  
  while (i <= length(Model)) {
    j <- 1
    
    # R2 / spearman correlation
    
    while (j <= length(Model[[i]])) {
      R2[["sites"]][[length(R2[["sites"]]) + 1]] <-
        names(Model[[i]])[j]
      
      R2[["locations"]][[length(R2[["locations"]]) + 1]] <-
        Model[[i]][[j]]$location
      
      R2[["people"]][[length(R2[["people"]]) + 1]] <-
        Model[[i]][[j]]$model$glmnet.fit$nobs
      
      R2[["valuesR2"]][[length(R2[["valuesR2"]]) + 1]] <-
        summary(lm(as.numeric(predict[[i]][[j]][, 2]) ~ as.numeric(Model[[i]][[j]][["expression"]][, 2])))$r.squared
      
      R2[["valuesSPM"]][[length(R2[["valuesSPM"]]) + 1]] <-
        cor(as.numeric(predict[[i]][[j]][, 2]), as.numeric(Model[[i]][[j]][["expression"]][, 2]), method = "spearman")
      
      
      # retrieve observed abd predicted values
      
      predicted[names(Model[[i]])[j]] <-
        list(as.numeric(predict[[i]][[j]][, 2]))
      
      observed[names(Model[[i]])[j]] <-
        list(as.numeric(Model[[i]][[j]][["expression"]][, 2]))
      
      
      
      
      
      j <- j + 1
      
    }
    
    i <- i + 1
    
  }
  
  if (type == 1) {
    # set bins and colour pallete for gradient depending on number of individuals
    
    bins <- cut(R2[["people"]], 10, include.lowest = TRUE)
    
    data <-
      data.frame(locations = R2[["locations"]],
                 values = R2[[paste("values", test, sep = "")]],
                 Data_Points = bins)
    
    p <- ggplot(data, aes(locations, values, color = Data_Points)) +
      geom_point() + scale_color_grey(start = 0.8, end = 0.2) + theme_classic()
    
    print(p)
    
  } else if (type == 2) {
    data <-
      data.frame(predicted = array(as.numeric(unlist(predicted))),
                 observed = array(as.numeric(unlist(observed))))
    
    p <- ggplot(data, aes(predicted, observed)) +
      geom_point()
    
    print(p)
    
    
  } else if (type == 3) {
    i <- 1
    
    while (i <= length(sModel)) {
      j <- 1
      
      while (j <= length(sModel[[i]])) {
        spredict[[i]][[j]] <-
          spredict[[i]][[j]][which(spredict[[i]][[j]][, 1] %in% sModel[[i]][[j]][["expression"]][, 1]), ]
        
        j <- j + 1
        
      }
      
      i <- i + 1
      
    }
    
    sR2 <-  list()
    
    i <- 1
    
    while (i <= length(sModel)) {
      j <- 1
      
      # R2 / spearman correlation
      
      while (j <= length(sModel[[i]])) {
        sR2[["sites"]][[length(sR2[["sites"]]) + 1]] <-
          names(sModel[[i]])[j]
        
        sR2[["locations"]][[length(sR2[["locations"]]) + 1]] <-
          sModel[[i]][[j]]$location
        
        sR2[["people"]][[length(sR2[["people"]]) + 1]] <-
          sModel[[i]][[j]]$model$glmnet.fit$nobs
        
        sR2[["valuesR2"]][[length(sR2[["valuesR2"]]) + 1]] <-
          summary(lm(as.numeric(spredict[[i]][[j]][, 2]) ~ as.numeric(sModel[[i]][[j]][["expression"]][, 2])))$r.squared
        
        sR2[["valuesSPM"]][[length(sR2[["valuesSPM"]]) + 1]] <-
          cor(as.numeric(spredict[[i]][[j]][, 2]),
              as.numeric(sModel[[i]][[j]][["expression"]][, 2]),
              method = "spearman")
        
        
        j <- j + 1
        
      }
      
      i <- i + 1
      
    }
    
    # get common points
    
    R2c <- list()
    sR2c <- list()
    
    R2c$valuesR2 <- R2$valuesR2[which(R2$sites %in% sR2$sites)]
    R2c$valuesSPM <- R2$valuesSPM[which(R2$sites %in% sR2$sites)]
    R2c$people <- R2$people[which(R2$sites %in% sR2$sites)]
    
    sR2c$valuesR2 <- sR2$valuesR2[which(sR2$sites %in% R2$sites)]
    sR2c$valuesSPM <- sR2$valuesSPM[which(sR2$sites %in% R2$sites)]
    
    bins <- cut(R2c[["people"]], 10, include.lowest = TRUE)
    
    data <-
      data.frame(
        Model_values = R2c[[paste("values", test, sep = "")]],
        sModel_values = sR2c[[paste("values", test, sep = "")]],
        Data_Points_Model = bins
      )
    
    p <- ggplot(data, aes(sModel_values, Model_values, color = Data_Points_Model)) +
      geom_point() + scale_color_grey(start = 0.8, end = 0.2) + theme_classic()
    
    print(p)
    
    
  }
  
  cat("Plot created, please check the plot window in RStudio")
  
}