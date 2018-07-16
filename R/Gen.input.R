#' Generate model input 
#'
#' This function allows you to output .RData files to be used as input for the analysis run in the second function.
#' @param Chromosome Specify chromosome 
#' @param ASE_file Specify file containing allele-specific expression information
#' @param legend_file Specify legend file containing SNP IDs and position info etc
#' @param haplotypes_file Provide haplotypes corresponding to SNPs in legend file, for given individuals found in samples file
#' @param Samples File containing cohort information
#' @param output_path Directory path to save output files to
#' @param Species Select same species as cohort. Defaults to "hsapiens" 
#' @param EnsemblVersion Specify version of Ensembl to download from. Defaults to NULL which is the most recent build 
#' @export
#' @examples
#' #Downloading Ensembl information for chromosome 22 of hsapiens, using the build corresponding to the sample data
#' Gen.input(22, "path_to_ASE_file/ASEfile.txt.gz",
#'                     "path_to_legend_file/file.legendfile.gz",
#'                    "path_to_haplotypes_file/hapfile.hap.gz",
#'                     "path_to_samples_file/samples.txt", EnsemblVersion=78)
#'
#' #Downloading Ensembl information for chromosome 22 of hsapiens, using the most recent build 
#' Gen.input(22, "path_to_ASE_file/ASEfile.txt.gz",
#'                     "path_to_legend_file/file.legendfile.gz",
#'                    "path_to_haplotypes_file/hapfile.hap.gz",
#'                     "path_to_samples_file/samples.txt")
#'
#' #Downloading Ensembl information for chromosome 1 of mmusculus, using the most recent build
#' Gen.input(1, "path_to_ASE_file/ASEfile.txt.gz",
#'                     "path_to_legend_file/file.legendfile.gz",
#'                   "path_to_haplotypes_file/hapfile.hap.gz",
#'                     "path_to_samples_file/samples.txt",
#'                     "mmusculus")



# 1. What you put into the command line
Gen.input <- function(Chromosome, ASE_file, legend_file, haplotypes_file, samples_file, output_path,
                      Species="hsapiens", EnsemblVersion=NULL)
{
  # 2. For timing length of script                                                                                                  
  Script.start.time <- Sys.time()
  # 3. Set your arguments. I will put a command that shows people what they have their arguments set to
  cat("\n\nARGUMENTS SET:\n")
  cat(paste(c("              Chromosome:      ", Chromosome, "\n"), collapse=""))
  cat(paste(c("              ASE file:        ", ASE_file, "\n"), collapse=""))
  cat(paste(c("              legend file:     ", legend_file, "\n"), collapse=""))
  cat(paste(c("              haplotypes file: ", haplotypes_file, "\n"), collapse=""))
  cat(paste(c("              Samples file:    ", samples_file, "\n"), collapse=""))
  #cat(paste(c("              Number of tasks: ", NumTasks, "\n"), collapse=""))
  cat(paste(c("              Species:         ", Species, "\n"), collapse=""))
  if(!is.null(EnsemblVersion))
  {
    cat(paste(c("              ensembl version:  ", EnsemblVersion, "\n"), collapse=""))
  } else{    
    cat(paste(c("              ensembl version: ", "latest build", "\n"), collapse=""))
  }    
  
  # 5. Set your results dataframe
  colHead<-c("id", "end", "Ind", "colnames(stat)", "(Intercept)_stat", "Reads_stat", "SEXmale_stat", "POPFIN_stat", "POPGBR_stat", "POPTSI_stat", "Variant_stat", "(Intercept)_p", "Reads_p", "SEXmale_
               p", "POPFIN_p", "POPGBR_p", "POPTSI_p", "Variant_p", "TSS", "Gene.x")
  
  # 6. Load input files                                 
  Samples<-readInputs(samples_file, "samples")
  LEG<-readInputs(legend_file, "legend")
  ASE<-readInputs(ASE_file, "allele counts")
  hap<-readInputs(haplotypes_file, "haplotypes")
  #ADD code to remove sites in ASE site not recorded as hets in hap file.
  
  if (nrow(LEG) != nrow(hap))
  {
    stop("Haplotype and legend files do not contain the same number of rows")
  }
  if ((ncol(hap)/2) != nrow(Samples))
  {
    stop("Must be two haplotypes per sample individual i.e. twice as many columns in the haplotype file as individuals in the samples file.")
  }
  colnames(hap)<-rep(Samples$ID,each=2)
  rownames(hap)<-LEG$id
  
  ######SHOULDNT ENFORCE THIS COLUMN
  LEG<-LEG[which(LEG$TYPE=="Biallelic_SNP"),]
  colnames(LEG)[c(2:4)]<-c("end", "ref", "alt")
  ASE_vars<-merge(ASE,LEG,by=c("end", "ref", "alt"), all.x=FALSE, all.y=FALSE)
  ASE_vars<-merge(ASE_vars,Samples, by.x="Ind", by.y="ID", all.x=FALSE, all.y=FALSE)
  #ASE_vars<-merge(ASE_vars, counts, by.x="Ind", by.y="ID", all.x=FALSE, all.y=FALSE)
  #some SNPs appear twice in ASE input when overlapping two genes
  ASE_vars<-ASE_vars[!duplicated(ASE_vars[,c("Ind", "id")]),]
  
  
  ensembl<-ensemblAttributes(Species, EnsemblVersion)
  
  
  
  # To correct for overall level of expression at given haplotype for given individual
  refSums<-cbind(aggregate(refCount~Ind, data=ASE, FUN=sum),0)
  altSums<-cbind(aggregate(altCount~Ind, data=ASE, FUN=sum),1)
  colnames(refSums)[2:3]<-c("Reads", "All")
  colnames(altSums)[2:3]<-c("Reads", "All")
  counts<-rbind(refSums,altSums)
  
  # MAY NEED TO CHANGE COLUMN NAMES. COULD JUST DO NUMBERS AND SAY WHAT EACH NUMBER CORRESPONDS TO IN A HELP FILE
  ensembl$TSS <- ifelse(ensembl$strand ==1, ensembl$transcript_start, ensembl$transcript_end)
  ensemblmerge <-   GRanges(seqnames= ensembl$chromosome_name,
                            ranges = IRanges(start=ensembl$exon_chrom_start, end=ensembl$exon_chrom_end, names=ensembl$ensembl_transcript_id),
                            geneid=ensembl$ensembl_gene_id, TSS=ensembl$transcript_start, strand= ensembl$strand)
  
  ASE_varsmerge <-  GRanges(seqnames=
                              Rle(ASE_vars$chr),
                            ranges = IRanges(start=ASE_vars$end, end=ASE_vars$end), names=ASE_vars$id)
  
  merge <- mergeByOverlaps(ASE_varsmerge,ensemblmerge)
  # MAY NEED TO CHANGE COLUMN NAMES. COULD JUST DO NUMBERS AND SAY WHAT EACH NUMBER CORRESPONDS TO IN A HELP FILE
  
  #MAYBE SHOULD MAKE DOING ASE SITES OUTSIDE A KNOWN GENE A SEPARATE OPTION
  #find out which snps are not found in my merge and so are without a known tss
  find_snps_without_tss <- which(!(ASE_vars$id %in% merge$names))
  #Pull them out of the metric file which also has their positions
  snps_without_tss <- ASE_vars[find_snps_without_tss,]
  #Do another merge in the same form as the previous one
  ensemblmerge <-   GRanges(seqnames= Rle(snps_without_tss$chr),
                            ranges = IRanges(start=snps_without_tss$end, end=snps_without_tss$end, names=snps_without_tss$id),
                            geneid=snps_without_tss$id, TSS=snps_without_tss$end, strand= snps_without_tss$strand)
  ASE_varsmerge <-  GRanges(seqnames=
                              Rle(snps_without_tss$chr),
                            ranges = IRanges(start=snps_without_tss$end, end=snps_without_tss$end), names=snps_without_tss$id)
  #Combine the snps without tss with those of a tss to test them all later
  merge2 <- mergeByOverlaps(ASE_varsmerge,ensemblmerge)
  completemerge <- rbind (merge,merge2)
  colnames(completemerge)[2] <- "id"
  # Add the TSS and geneid information to the ASE_vars df
  newdf <- data.frame(id=completemerge$id, TSS=completemerge$TSS, Gene=completemerge$geneid)
  # This merge takes a while
  cat("Now merging\n")
  ASE_vars <- merge(newdf, ASE_vars, by = "id")
  ASE_vars<-ASE_vars[!duplicated(ASE_vars[,c("Ind", "id")]),]
  ASE_vars$propRef<-ASE_vars$refCount/(ASE_vars$refCount+ASE_vars$altCount)
  ASE_vars$totalReads<-rowSums(ASE_vars[,c("refCount","altCount")])
  ASE_vars$logRatio<-log2((ASE_vars$refCount+1)/(ASE_vars$altCount+1))
  
  expectedRatio <- 0.5
  if(length(ASE_vars$propRef > 1000))
  {
    expectedRatio <- median(ASE_vars$propRef)
  }
  cat("Setting expected proportion of reference reads to ", expectedRatio, "in binomial test\n")
  ASE_vars$binomp<-mapply(applyBinom, ASE_vars$refCount, ASE_vars$refCount+ASE_vars$altCount, expectedRatio)
  cat("Merge complete\n")
  output_file = paste(c(output_path, "Run.model.input_Chr",  Chromosome, ".RData"), collapse="")
  cat(paste(c("Saving ", output_file, "\n"), collapse=""))
  
  save(list = ls(all.names = TRUE), file = output_file, envir = environment())
  cat("Finished\n")
}  


applyBinom<-function (x,n, p){
  binom.test(x,n, p)$p.value
}


readInputs<-function(thisFile, type)
{
  if(file.exists(thisFile))
  {
    if(type == "haplotypes")
    {
      thisData<-read.table(thisFile, header=F)  
    } else {
      thisData<-read.table(thisFile, stringsAsFactors = F, header=T)
    }
    cols <- colnames(thisData)
    #lets change samples file so can take any number of different covariates
    #do we need the ref and alt alleles in ASE file? Best to keep required columns to a minimum
    #also lets call allele counts rather than ASE
    if((type=="samples") & (cols[1] != "ID")) 
    {
      print(sample_example)
      stop("First line of samples file does not start with ID. Make sure a header line is specified as above")
      
    } else if((type=="allele counts") & !all(colnames(ASE_example) %in% cols)) {
      print(ASE_example)
      stop("input file does not contain the correct column names. Allele counts file must be in the format shown above.")
      
    } else if((type=="legend") & !all(colnames(LEG_example) %in% cols)) {
      
      print(LEG_example)
      stop("input file does not contain the correct column names. Legend file must be in the format shown above.")
      
    } else if(type=="haplotypes" & (length(which(thisData != 0 & thisData != 1)) > 0)) {
      stop("Error: Haplotype file contains values other than 0 or 1.")
    } else {
      cat(type, "file successfully loaded\n")
    }
  } else {
    stop(paste("Couldnt find specified ", type," file at: ", getwd(), "/", thisFile, sep=""))
  }
  return(thisData)
}

