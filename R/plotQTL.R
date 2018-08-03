plotQTL<-function(modelOut, aseVar, qtlVar)
{
  qtlData<-mergeExprGenos(aseVar, modelOut, qtlVar)
  exprData<-qtlData$exprVar
  exprData$RPM<-exprData$Count/(exprData$Reads/1000000)
 
  dodge <- position_dodge(width = 0.7)
  p1<-ggplot(exprData, aes(x=as.factor(get(qtlVar)), y=RPM))+
    geom_violin(trim=FALSE, position = dodge, aes(fill = as.factor(get(qtlVar))))+
    geom_boxplot(width=0.2, position = dodge, fill = "white")+ 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
    ylab(paste("RPM at ", aseVar, sep="")) + xlab(paste("Allele at ", qtlVar,sep="")) +
    theme_pubr() + stat_compare_means(label.x = 1.3)+ theme(legend.position="none")
  
  p1
} 