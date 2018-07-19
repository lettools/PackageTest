plotASEMetrics<-function(input, individuals=NULL, genes=NULL, variants=NULL)
{
  thisASE<-input$ASE
  if(!is.null(individuals))
  {
    thisASE<-input$ASE[which(input$ASE$Ind %in% individuals),]
  }
  if(!is.null(genes))
  {
    thisASE<-input$ASE[which(input$ASE$Gene.y %in% genes),]
  }
  if(!is.null(variants))
  {
    thisASE<-input$ASE[which(input$ASE$ID %in% variants),]
  }
  if(dim(input$ASE)[1] > 0)
  {
    numInd<-length(unique(thisASE$Ind))
    numVar<-length(unique(thisASE$ID))
    numGenes<-length(unique(thisASE$Gene.y[!is.na(thisASE$Gene.y)]))
    numNA<-length(unique(thisASE$ID[is.na(thisASE$Gene.y)]))
    cat("Numbers after filtering:\n")
    cat("\t", numInd, "individuals\n")
    cat("\t", numVar, "variants\n")
    cat("\t", numGenes, "genes\n")
    
    title<-paste("ASE metrics across ", numInd, " individuals, ", numGenes, " genes and ", numVar, " variants (", numNA, " outside known genes)", sep="")
    
    hetCounts<-sort(table(thisASE$ID), decreasing = TRUE)
    #check plotting when only one SNP
    hetCounts<-as.vector(hetCounts[which(hetCounts>0)])
    hetCountsRank<-cbind.data.frame(Individuals=hetCounts, Rank=1:length(hetCounts))
    colnames(hetCountsRank)<-c("Individuals.Freq", "Rank")
    p1<-ggplot(hetCountsRank, aes(x=Rank, y=Individuals.Freq)) +
      geom_point() + scale_y_continuous(trans='log10') + annotation_logticks(scaled = TRUE,sides = "lr") + 
      xlab("Rank of variant") + ylab("Number of heterozygote individuals") + theme_pubr()
    
    #rather convoluted command to get variant counts by gene
    geneCounts<-table(as.character(unique(thisASE[complete.cases(thisASE[,c("ID","Gene.y")]),c("ID","Gene.y")])$Gene.y))
    geneCountsRank<-cbind.data.frame(Genes=sort(geneCounts, decreasing = TRUE), Rank=1:length(geneCounts))
    p2<-ggplot(geneCountsRank, aes(x=Rank, y=Genes.Freq)) +
      geom_point() + scale_y_continuous(trans='log10') + annotation_logticks(scaled = TRUE,sides = "lr") + 
      xlab("Rank of gene") + ylab("Number of heterozygote variants") + theme_pubr()
    
    p3<-ggplot(thisASE, aes(x=propRef)) +
      geom_histogram(aes(y=..density..),binwidth=.05, colour="black", fill="white") +
      geom_vline(aes(xintercept=median(propRef, na.rm=T)),   # Ignore NA values for median
                 color="red", linetype="dashed", size=1) + xlab("Proportion of reads carrying reference allele") +
      geom_density(alpha=.2, fill="#FF6666") + theme_pubr()
    
    logTransBi<--log10(thisASE$binomp)
    maxLogP<-max(logTransBi[is.finite(logTransBi)])
    p4<-ggplot(thisASE, aes(x=logRatio, y=totalReads)) +
      geom_point(aes(colour=-log10(binomp+1e-300))) + theme_pubr() + 
      xlab("Log2 ratio ((Ref. reads + 1)/(Alt. reads +1))") + ylab("Total number of reads") +
      scale_y_continuous(trans='log10') +
      scale_colour_gradientn(name = "-log10(Binomial P value)",colors=c("cornflowerblue","orange", "red"), values=c(0,3/maxLogP,1))
      
    annotate_figure(ggarrange(p1,p2,p3,p4), top=title)
  } else {
    stop("No data left to plot. Did you specify valid ids?")
  }
  
}