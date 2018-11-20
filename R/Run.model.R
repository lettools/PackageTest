#' Run Model
#'
#' Allows you to measure statistical association between nearby regulatory variants and the level of expression at a heterozygous coding polymorphism, controlling for factors such as sex and population,
#' by utilising a generalized linear model and applying permutations to the data in order to provide a robust p-value
#' @param input_file This is the RData file outputted from the first function, gen.input
#' @param Task This analysis is very computationally burdensome. To speed up the process it is an advantage to split it up into tasks that may be run on multiple nodes, concurrently. 
#' @param progress_path The function allows checkpointing, which allows the program to be killed and picked up again at a later stage, without starting from the beginning. Occurs every 2 hours
#' @param numTasks Select how many jobs/tasks to split the process up into. Defaults to 100
#' @param Chromosome Specify chromosome used in input file
#' @param numPerms Select how many permutations to run. Along with splitting the process up into simultaneous tasks, this is the biggest factor in determining how long the analysis will take. However, the 
#'more permutations, in general and up until a point, the more precise and accurate the results may be; for example, if set to 100, the minimum p-value that can possibly be reached as a result of 
#'permutations, is 0.01. Defaults to 100,000. 
#' @param TSSwindow This represents the distance over which nearby variants will be selected,  either side of the transcript start site. Defaults to 500kb 
#' @param pval_threshold There is a theoretical minimum p-value for each particular combination of reference and alternative alleles for a given set of individuals for a given nearby variant of an ASE site
#'. 
#' \t\t\t This parameter sets the upper limit. Default is 0.00005. In this example the model will not be run if it is not possible to reach a p-value as low as 0.00005, even theoretically. 
#' @export
#' @examples
#' #Run model with task set to 10, chromosome to 22, for 100,000 permutations, a transcript start site window of 500kb and a theoretical p-value threshold of 0.00005
#' Run.Model('input_file.RData', 10, 
#'         'path_to_progress_file',
#'         'output_path',
#'         Chromosome=22)
#' 
#' #Run model with task\tset to 10, chromosome to 22, for 10,000 permutations, a transcript start site window of 500kb and a theoretical p-value of 0.0005
#' Run.Model('input_file.RData', 10,
#'         'path_to_progress_file',
#'         'output_path',
#'         Chromosome=22, numPerms=10000, pval_threshold=10000)
#'
#' #Run model with task\tset to 2, chromosome to 12, for 10,000 permutations, a transcript start site window of 1Mb and a theoretical p-value of 0.0005
#' Run.Model('input_file.RData', 2,
#'         'path_to_progress_file',
#'         'output_path',
#'         Chromosome=22, numPerms=10000,TSSwindow=100000, pval_threshold=10000)
#'
#'NB: The smallest possible p-value attainable as a result of running permutations is 1/numPerms. Hence, there is no advantage to setting the minimum p-value threshold to below this number.


Run.Model <- function(inputObj, output_prefix, task = 1, totalTasks = 1, minInd = 10, numPerms = 1e+05, TSSwindow = 5e+05, pval_threshold = 5e-05, other_all = FALSE) {
    # progress_file = paste(c(progress_path, 'progress_', Chromosome, '_task', Task, '.RData'), collapse='')
    
    
    if (is.numeric(task) & is.numeric(totalTasks) & (task <= totalTasks)) {
        hetCounts <- getHetCounts(task, totalTasks, minInd, inputObj)
        
    } else {
        stop("Task or total tasks numbers are misspecified")
    }
    
    
    
    
    
    # i is the ASE SNP currently on. This changes if a progress file has a different value for it
    i <- 1
    results <- data.frame()
    task.start.time <- Sys.time()
    cycle.start.time <- Sys.time()
    
    # This is a checkpoint of results and which current ASE SNP on.
    progress_file = paste(c(output_prefix, "_progress_", task, "_", totalTasks, ".RData"), collapse = "")
    # progress_file <- paste(c(progress_path, 'progress_Chr', Chromosome, '_task', task, '.RData'), collapse = '') Check if the
    # progress_file exists.  If it does, check that it is above a size of zero, or it will produce errors If it doesn't, create it.
    if (file.exists(progress_file)) {
        if (file.size(progress_file) > 0) {
            load(file = progress_file)
        }
    } else {
        file.create(progress_file)
    }
    # cat(paste(c('At the beginning of this cycle, the results df stored in progress_file now has ', dim(results)[1], ' lines.\n'),
    # collapse = ''))
    
    
    covarNames <- colnames(inputObj$covar)
    
    # Start of model loop
    if (dim(hetCounts)[1] > 0) {
        while (i <= dim(hetCounts)[1]) {
            # i <- i + 1
            cat(paste(c("Currently on ASE SNP", i, "of", dim(hetCounts)[1], "\n"), collapse = " "))
            i.start.time <- Sys.time()
            cycle.end.time <- Sys.time()
            cycle.time.diff <- difftime(cycle.end.time, cycle.start.time, units = "hours")
            if (cycle.time.diff[1] >= 2) {
                # Reset the cycle start time
                cycle.start.time <- Sys.time()
                # Save where you are in hetCounts, so you know which ASE SNPs you've done so far in this vector, and save the results df
                save(i, results, file = progress_file)
                cat("results df saved to progress_file.\n")
            }
            
            # need to think about these two lines. Does it mean same variant will be run multiple times?
            thisId <- hetCounts$ID[i]
            nearbyVar <- inputObj$leg$ID[which(inputObj$leg$end > (hetCounts$TSS[i] - TSSwindow) & inputObj$leg$end < (hetCounts$TSS[i] + 
                TSSwindow))]
            mergeResults <- mergeExprGenos(thisId, inputObj, nearbyVar)
            
            exprVar <- mergeResults$exprVar
            beforenumber_of_nearby_variants <- dim(exprVar)[2]
            if (dim(exprVar)[[1]] < 171) {
                failCols <- mergeResults$min_pval[which(mergeResults$min_pval > pval_threshold)]
                exprVar <- exprVar[, which(!(colnames(exprVar) %in% names(failCols)))]
            }
            # Changed this to implement min-pvalue threshold
            number_of_nearby_variants <- dim(exprVar)[2]
            cat(paste(c("This coding variant has", number_of_nearby_variants, "of", beforenumber_of_nearby_variants, "that can theoretically pass the minimum permutation p-value threshold.\n"), 
                collapse = " "))
            number_individuals <- dim(exprVar)[1]/2
            cat(paste(c("This coding variant is heterozygote at", number_individuals, "individuals.\n"), collapse = " "))
            
            theseResults <- fitModels(exprVar, covarNames, hetCounts[i, ],other_all)
            # print('theseResults') print(head(theseResults))
            if (dim(theseResults)[1] > 0) {
                cat("theseResults df exists.\n")
                theseResults$numPerm <- 0
                theseResults$numPermExceed <- 0
                
                if (numPerms > 0) {
                  # do the permutations
                  perms <- 100
                  totalPerms <- 100
                  cat("Beginning permutations\n")
                  # print(Sys.time())
                  while (totalPerms <= numPerms) {
                    # get just nominal signif
                    toPerm <- theseResults[which(theseResults$Variant_p < 5e-06 & theseResults$numPermExceed < 5), ]
                    numLeft <- dim(toPerm)[1]
                    # cat(paste(c(totalPerms, '\n'), collapse=' '))
                    
                    if (numLeft > 0) {
                      # set the variant coeff to numeric as should be no NAs now
                      # toPerm$Variant_coeff<-as.numeric(levels(toPerm$Variant_coeff))[toPerm$Variant_coeff]
                      exprVar2 <- exprVar[, c(covarNames, "All", "Count", "Reads", as.character(toPerm[, "colnames(stat)"]))]
                      permResults <- fitPerms(exprVar2, toPerm, perms, covarNames, hetCounts[i, ],other_all)
                      theseResults[rownames(theseResults) %in% rownames(permResults), ]$numPerm <- permResults$numPerm
                      theseResults[rownames(theseResults) %in% rownames(permResults), ]$numPermExceed <- permResults$numPermExceed
                    }
                    print(paste(c(totalPerms, " permutations completed"), collapse = ""))
                    totalPerms <- totalPerms + perms
                  }
                }
                results <- rbind.fill(results, theseResults[which(theseResults$Variant_p <= 1), ])
            }
            # if increment here, only increments if the loop completes.
            i.end.time <- Sys.time()
            i.time.diff <- difftime(i.end.time, i.start.time, units = "hours")
            cat(paste(c("ASE SNP took", i.time.diff[[1]], "hours to complete\n"), collapse = " "))
            i <- i + 1
        }
    } else {
        stop("No coding heterozygote variants to analyse in this task\n")
    }
    
    output_file = paste(c(output_prefix, "_", task, "_", totalTasks, ".txt"), collapse = "")
    write.table(results, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
    
    inputObj$pvalues <- results
    cat("Task finished\n")
    output_file = paste(c(output_prefix, "_", task, "_", totalTasks, ".RData"), collapse = "")
    cat(paste(c("Saving ", output_file, "\n"), collapse = ""))
    save(inputObj, file = output_file)
    return(inputObj)
    # Script.end.time <- Sys.time() Script.duration <- difftime(Script.end.time, Script.start.time, units='hours') cat(paste(c('Script
    # took ', Script.duration[1], ' minutes\n'), collapse=''))
}

getHetCounts <- function(set, setTotal, numInd, inputO) {
    # aggregate ignores rows with NAs so use dplyr hetCountsAll <- aggregate(Ind ~ ID + end + TSS + Gene, data = inputO$ASE, FUN = length)
    hetCountsAll <- group_by(inputO$ASE, ID, end, TSS, Gene) %>% summarise(count = length(Ind))
    hetCountsPass <- hetCountsAll[which(hetCountsAll$count >= numInd), ]
    eachTSS <- unique(hetCountsPass$TSS)
    if (set > length(eachTSS)) {
        stop("Task number is greater than number of distinct TSSs")
    }
    toKeep <- eachTSS[seq(set, length(eachTSS), setTotal)]
    subHetCounts <- hetCountsPass[which(hetCountsPass$TSS %in% toKeep), ]
    cat(paste(c("A total of", dim(hetCountsAll)[1], "coding heterozygote variants found", "\n"), collapse = " "))
    cat(paste(c(dim(hetCountsPass)[1], "coding heterozygote variants found in at least", numInd, "individuals\n"), collapse = " "))
    cat(paste(c("Analysing", dim(subHetCounts)[1], "coding heterozygote variants in this task", set, "of", setTotal, "\n"), collapse = " "))
    return(subHetCounts)
}

# Set functions
getFitStats <- function(fits, theseHetCounts) {
    coeff_p <- data.frame()
    # remove model fits that failed then check some are left
    fits_NN <- fits[!sapply(fits, is.null)]
    # variant alleles can be perfectly correlated with another variable leading to NAs. Exclude these here
    fits_NN <- fits_NN[!sapply(fits_NN, function(f) {
        a <- f$coefficients
        is.na(a[length(a)])
    })]
    
    if (length(fits_NN) > 0) {
        stat <- sapply(fits_NN, function(f) {
            c <- summary(f)$coefficients[, 3]
            if (is.factor(c)) {
                as.numeric(levels(c))[c]
            } else {
                c
            }
        })
        pvals <- sapply(fits_NN, function(f) {
            p <- summary(f)$coefficients[, 4]
            if (is.factor(p)) {
                as.numeric(levels(p))[p]
            } else {
                p
            }
        })
        # get names of variants that didnt return null
        testedVar <- sapply(fits_NN, function(f) {
            attributes(summary(f)$terms)$variables[[length(attributes(summary(f)$terms)$variables)]]
        })
        testedVar<-gsub("^sample\\(|\\)$", "", testedVar)
        colnames(stat) <- testedVar
        colnames(pvals) <- testedVar
        #tidy up the row names
        rownames(stat) <- paste(names(coef(fits_NN[[1]])), "_stat", sep = "")
        rownames(pvals) <- paste(names(coef(fits_NN[[1]])), "_p", sep = "")
        rownames(stat) <- gsub("sample\\(|\\)", "", row.names(stat))
        rownames(pvals) <- gsub("sample\\(|\\)", "", row.names(pvals))
        rownames(stat) <- gsub(testedVar, "Variant", row.names(stat))
        rownames(pvals) <- gsub(testedVar, "Variant", row.names(pvals))
        rownames(stat) <- gsub("\\[seq_len\\(length\\(Variant\\)\\) \\+ c\\(1, -1\\)\\]", "_alt_all", row.names(stat))
        rownames(pvals) <- gsub("\\[seq_len\\(length\\(Variant\\)\\) \\+ c\\(1, -1\\)\\]", "_alt_all", row.names(pvals))
        
        #grepping for seq_len may cause problems in the admittedly rare occurrence of a user having a covariate with this in
        #so should look at updating this.
        #rownames(stat)[grepl( "seq_len" , rownames( stat ) )]<-"Other_allele_stat"
        #rownames(pvals)[grepl( "seq_len" , rownames( pvals ) )]<-"Other_allele_p"
        #rownames(stat)[length(rownames(stat))] <- "Variant_stat"
        #rownames(pvals)[length(rownames(pvals))] <- "Variant_p"
        coeff_p <- cbind.data.frame(theseHetCounts, cbind.data.frame(colnames(stat), cbind.data.frame(t(stat), t(pvals))))
        # add NAs for any missing columns NEED TO FIX WHAT HAPPENS IF COLUMNS MISSING missingCols <- colHead[which(!(colHead %in%
        # colnames(coeff_p)))] coeff_p[, missingCols] <- 'NA' order columns coeff_p <- coeff_p[, colHead]
    }
    return(coeff_p)
}

fitPerms <- function(exprGenos, nomResults, numPerm, covarNames, theseHetCounts, altAll) {
    
    # order nomResults so that can reorder after merge to maintain row positions COMMENTED OUT NEXT LINE
    # nomResults<-nomResults[order(nomResults$'colnames(coeff)'),]
    for (k in 1:numPerm) {
        vars <- colnames(exprGenos)[-c(1:(length(covarNames) + 3))]
        
        if(isTRUE(altAll) {
          fitsA <- lapply(vars, function(x) {
            frm <- as.formula(paste(paste("Count ~ Reads", paste(covarNames[-1], collapse = " + "),paste(substitute(j, list(j = as.name(x))),"[seq_len(length(",substitute(j, list(j = as.name(x))),")) + c(1,-1)] * ", sep=""), " sample(", sep = " + "), substitute(k, list(k = as.name(x))), ")", sep = ""))
            tryCatch(glm.nb(frm, data = exprGenos), error = function(e) NULL)
          })
        } else {
        
          fitsA <- lapply(vars, function(x) {
              frm <- as.formula(paste(paste("Count ~ Reads", paste(covarNames[-1], collapse = " + "), "sample(", sep = " + "), substitute(k, 
                  list(k = as.name(x))), ")", sep = ""))
              tryCatch(glm.nb(frm, data = exprGenos), error = function(e) NULL)
          })
        }
        fitStats <- getFitStats(fitsA, theseHetCounts)
        if (dim(fitStats)[1] > 0) {
            fitStats$"colnames(stat)" <- gsub("sample\\(|\\)", "", fitStats$"colnames(stat)")
            
            # get rows where have a valid permutation count so can increment just their counts
            a <- which(nomResults$"colnames(stat)" %in% fitStats$"colnames(stat)")
            nomResults$numPerm[a] <- nomResults$numPerm[a] + 1
            
            
            mergedStat <- join(nomResults[, c("colnames(stat)", "Variant_stat")], fitStats[, c("colnames(stat)", "Variant_stat")], by = "colnames(stat)", 
                type = "left")
            
            # get rows where the permutation coefficent is greater than the nominal coefficent so can increment their counts
            nomResults$numPermExceed[which(abs(mergedStat[, 3]) >= abs(mergedStat[, 2]))] <- nomResults$numPermExceed[which(abs(mergedStat[, 
                3]) >= abs(mergedStat[, 2]))] + 1
            
            
        }
    }
    return(nomResults)
}

fitModels <- function(exprGenos, covarNames, theseHetCounts, altAll) {
    vars <- colnames(exprGenos)[-c(1:(length(covarNames) + 3))]
    
    #frm <- as.formula(paste("Count ~ Reads", paste(covarNames[-1], collapse = " + "), sep=" + "))
    #residFit<-glm.nb(frm, data=exprGenos)
    #exprGenos$residuals<-exprGenos$Count - predict(residFit)
    
    #if should fit allele on other chromosome as well
    if(isTRUE(altAll)) {
      fitsA <- lapply(vars, function(x) {
        frm <- as.formula(paste("Count ~ Reads", paste(covarNames[-1], collapse = " + "), paste(substitute(j, list(j = as.name(x))),"[seq_len(length(",substitute(j, list(j = as.name(x))),")) + c(1,-1)] * ", substitute(j, list(j = as.name(x))), sep=""), sep = " + "))
        tryCatch(glm.nb(frm, data = exprGenos), error = function(e) NULL)
      })
    } else {
      fitsA <- lapply(vars, function(x) {
          frm <- as.formula(paste("Count ~ Reads", paste(covarNames[-1], collapse = " + "), substitute(j, list(j = as.name(x))), sep = " + "))
          tryCatch(glm.nb(frm, data = exprGenos), error = function(e) NULL)
      })
    }
    return(getFitStats(fitsA, theseHetCounts))
}

#test<-exprGenos %>% group_by(.dots=c("Ind", "POP", "SEX")) %>% summarise_all(funs(sum))
#fitsA[[1]]<-glm.nb(Count ~ Reads + POP + SEX + rs2017689, exprGenos)
#getFitStats(fitsA, theseHetCounts) 
  