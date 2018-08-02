plotQTL<-function(exprData, var)
{
  exprData$RPM<-exprData$Count/(exprData$Reads/1000000)
 
  dodge <- position_dodge(width = 0.7)
  p1<-ggplot(exprData, aes(x=as.factor(rs117966628), y=RPM))+
    geom_violin(trim=FALSE, position = dodge, aes(fill = as.factor(rs117966628)))+
    geom_boxplot(width=0.2, position = dodge, fill = "white")+ 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    ylab("RPM") + xlab("Allele") +
    theme_pubr() + stat_compare_means(label.x = 1.3)
  
  p1
} 