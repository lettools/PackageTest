#' Generate plots of expression values for given eQTL and expressed heterozygote variant
#'
#' This function takes the output of Run.model and produces various plots for the specified putative regulatory variant.
#' @param modelOut The output of the Run.model function
#' @param aseVar The id of the allele specific expression variant.
#' @param qtlVar The id of the putative regulatory variant.
#' @export
#' @examples
#' #' #Plot the output of Run.model.R stored in modelOut
#' plotQTL(aseDat,'rs116197803', 'rs175183')


plotQTL <- function(modelOut, aseVar, qtlVar) {
    qtlData <- mergeExprGenos(aseVar, modelOut, qtlVar)
    exprData <- qtlData$exprVar
    exprData$RPM <- exprData$Count/(exprData$Reads/1e+06)
    
    dodge <- position_dodge(width = 0.7)
    p1 <- ggplot(exprData, aes(x = as.factor(get(qtlVar)), y = RPM)) + geom_violin(trim = FALSE, position = dodge, aes(fill = as.factor(get(qtlVar)))) + 
        geom_boxplot(width = 0.2, position = dodge, fill = "white") + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) + 
        ylab(paste("RPM at ", aseVar, sep = "")) + xlab(paste("Allele at ", qtlVar, sep = "")) + theme_pubr() + stat_compare_means(label.x = 1.3) + 
        theme(legend.position = "none")
    
    p1
}
