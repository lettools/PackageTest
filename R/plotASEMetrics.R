plotASEMetrics<-function(individuals=NULL, genes=NULL, variants=NULL)
{
  thisASE<-ASE_vars
  if(!is.null(individuals))
  {
    thisASE<-ASE_vars[which(ASE_vars$Ind %in% individuals),]
  }
  if(!is.null(genes))
  {
    thisASE<-ASE_vars[which(ASE_vars$Gene.y %in% genes),]
  }
  if(!is.null(variants))
  {
    thisASE<-ASE_vars[which(ASE_vars$id %in% variants),]
  }
  if(dim(ASE_vars)[1] > 0)
  {
    numInd<-length(unique(thisASE$Ind))
    numVar<-length(unique(thisASE$id))
    numGenes<-length(unique(thisASE$Gene.y[!is.na(thisASE$Gene.y)]))
    numNA<-length(unique(thisASE$id[is.na(thisASE$Gene.y)]))
    cat("Numbers after filtering:\n")
    cat("\t", numInd, "individuals\n")
    cat("\t", numVar, "variants\n")
    cat("\t", numGenes, "genes\n")
    
    title<-paste("Allelic imbalances across ", numInd, " individuals, ", numGenes, " genes and ", numVar, " variants (", numNA, " outside known genes)", sep="")
    p1<-ggplot(thisASE, aes(x=propRef)) +
      geom_histogram(aes(y=..density..),binwidth=.05, colour="black", fill="white") +
      geom_vline(aes(xintercept=median(propRef, na.rm=T)),   # Ignore NA values for median
                 color="red", linetype="dashed", size=1) + xlab("Proportion of reads carrying reference allele") +
      geom_density(alpha=.2, fill="#FF6666") + theme_pubr()
    p2<-ggplot(thisASE, aes(x=logRatio, y=totalReads)) +
      geom_point() + theme_pubr() + xlab("Log2 ratio ((Ref. reads + 1)/(Alt. reads +1))") + ylab("Total number of reads") +
      scale_y_continuous(trans='log10')
    annotate_figure(ggarrange(p1,p2), top=title)
  } else {
    stop("No data left to plot. Did you specify valid ids?")
  }
  
}