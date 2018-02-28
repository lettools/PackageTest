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
Gen.input <- function(Chromosome, ASE_file, legend_file, haplotypes_file, Samples_file, output_path,
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
  cat(paste(c("              Samples file:    ", Samples_file, "\n"), collapse=""))
  #cat(paste(c("              Number of tasks: ", NumTasks, "\n"), collapse=""))
  cat(paste(c("              Species:         ", Species, "\n"), collapse=""))
  if(!is.null(EnsemblVersion))
  {
    cat(paste(c("              ensembl version:  ", EnsemblVersion, "\n"), collapse=""))
  } else{    cat(paste(c("              ensembl version: ", "latest build", "\n"), collapse=""))
}    
    
    # 4. Loading packages
    library(GenomicRanges)
    library(reshape2)
   # library(rtracklayer)
#    library(speedglm)
    library(stringr)
    library(plyr)
    library(biomaRt)
    cat("\nAll packages loaded\n")
    
    # 5. Set your results dataframe
    colHead<-c("id", "end", "Ind", "colnames(stat)", "(Intercept)_stat", "Reads_stat", "SEXmale_stat", "POPFIN_stat", "POPGBR_stat", "POPTSI_stat", "Variant_stat", "(Intercept)_p", "Reads_p", "SEXmale_
               p", "POPFIN_p", "POPGBR_p", "POPTSI_p", "Variant_p", "TSS", "Gene.x")
    
    # 6. Load input files                                 
    Samples<-read.table(Samples_file, header=T)
    cols <- colnames(Samples)
    if (!"ID" %in% cols | !"SEX" %in% cols | !"POP" %in% cols)
    {
      ID <- c("HG00096", "HG00097", "HG00099", "NA20827") 
      SEX <- c("male", "female", "female", "male")
      POP <- c("GBR", "GBR", "GBR", "TSI")
      sample_example <- data.frame(ID, SEX, POP) 
      print(sample_example)
      cat("\n")
      stop("input file does not contain the correct column names. Samples file must be in the format given above.")
    }
    print("Samples loaded")
    ASE<-read.table(ASE_file, stringsAsFactors = F, header=T)
    cols <- colnames(ASE)
    if (!"chr" %in% cols | !"start" %in% cols | !"end" %in% cols | !"ref" %in% cols | !"Ind" %in% cols | !"refCount" %in% cols | !"altCount" %in% cols)
    {
      chr <- c(22,22,22)
      start <- c(135031, 135031, 135031)
      end <- c(135032, 135032, 135032)
      ref <- c("G","G","G")
      alt <- c("A","A","A")
      Ind <- c("HG00276", "HG00282", "NA11831")
      refCount <- c(19,12,10)
      altCount <- c(0,0,0)
      ASE_example <- data.frame(chr, start, end, ref, alt, Ind, refCount, altCount)
      print(ASE_example)
      cat("\n")
      stop("input file does not contain the correct column names. ASE file must be in the format given above.")
    }            
    print("ASE file loaded")
    #legend_file <- "/exports/eddie/scratch/s1463821/1000GP_Phase3_chr22.EUR.noFixed.legend.gz" 
    LEG<-read.table(legend_file, header=T)
    cols <- colnames(LEG)
    if (!"id" %in% cols | !"position" %in% cols | !"a0" %in% cols | !"a1" %in% cols | !"TYPE" %in% cols)
    {
      id <- c("10:60515:C:T", "rs148087467:60523:T:G", "rs147855157:61372:CA:C")
      position <- c(60515, 60523, 61372)
      a0 <- c("C", "T", "CA")
      a1 <- c("T", "G", "C")
      TYPE <- c("Biallelic_SNP", "Biallelic_SNP", "Biallelic_INDEL")
      LEG_example <- data.frame(id, position, a0, a1, TYPE)
      print(LEG_example)
      cat("\n")
      stop("input file does not contain the correct column names. LEG file must be in the format given above.")
    }
    print("Legend file loaded")
    LEG$hapLine<-seq(1:dim(LEG)[1])
    # haplotypes_file <- "/exports/eddie/scratch/s1463821/1000GP_Phase3.EUR.chr22.noFixed.hap.gz"  
    hap<-read.table(haplotypes_file, header=F)
    print("haplotypes file loaded")
    if (nrow(LEG) != nrow(hap))
    {
      stop("Haplotype and legend files do not contain the same number of rows")
    }
    if ((ncol(hap)/2) != nrow(Samples))
    {
      stop("Must be two haplotypes per sample individual i.e. twice as many hap columns as Sample rows")
    }
    colnames(hap)<-rep(Samples$ID,each=2)
    rownames(hap)<-LEG$id
    LEG<-LEG[which(LEG$TYPE=="Biallelic_SNP"),]
    colnames(LEG)[c(2:4)]<-c("end", "ref", "alt")
    ASE_vars<-merge(ASE,LEG,by=c("end", "ref", "alt"), all.x=FALSE, all.y=FALSE)
    ASE_vars<-merge(ASE_vars,Samples, by.x="Ind", by.y="ID", all.x=FALSE, all.y=FALSE)
    #ASE_vars<-merge(ASE_vars, counts, by.x="Ind", by.y="ID", all.x=FALSE, all.y=FALSE)
    #some SNPs appear twice in ASE input when overlapping two genes
    ASE_vars<-ASE_vars[!duplicated(ASE_vars[,c("Ind", "id")]),]
    
    ensemblConnect<-function(species="hsapiens", ensemblVersion=NULL) {
      if(!is.null(ensemblVersion))
      {
        ensembl = useEnsembl(biomart="ensembl", dataset=paste(species, "_gene_ensembl", sep=""), version=ensemblVersion, host="uswest.ensembl.org")
      }
      else
      {
        ensembl = useEnsembl(biomart="ensembl", dataset=paste(species, "_gene_ensembl", sep=""), host="uswest.ensembl.org")
      }
      return(ensembl)
    }
    
    ens <- NULL
    attempt <- 0
    while( is.null(ens) && attempt <= 19 ) {
      attempt <- attempt + 1
      print(paste("Trying to connect to Ensembl. Attempt", attempt, sep=" "))
      tryCatch({
        # Obtain dataframe with ensembl information
        ens<-ensemblConnect()
      } , error = function(err) {
        print("Connection failed.")
      } , finally = {
        if(!is.null(ens))
        {
          print("Connection successful.")
        }
      }
      )
    } 
    ensembl <- getBM(attributes=c("chromosome_name", "ensembl_transcript_id", "ensembl_gene_id", "start_position", "exon_chrom_start", "exon_chrom_end", "transcript_start" , "transcript_end", "strand"), mart
                     =ens)
    print("Downloading ensembl data") 
    print("ensembl data downloaded")
    #print(head(ensembl))
    
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
                                Rle(rep(as.character(Chromosome), each=dim(ASE_vars)[1])),
                              ranges = IRanges(start=ASE_vars$end, end=ASE_vars$end), names=ASE_vars$id)
    
    merge <- mergeByOverlaps(ASE_varsmerge,ensemblmerge)
    # MAY NEED TO CHANGE COLUMN NAMES. COULD JUST DO NUMBERS AND SAY WHAT EACH NUMBER CORRESPONDS TO IN A HELP FILE
    
    #find out which snps are not found in my merge and so are without a tss
    find_snps_without_tss <- which(!(ASE_vars$id %in% merge$names))
    #Pull them out of the metric file which also has their positions
    snps_without_tss <- ASE_vars[find_snps_without_tss,]
    #Do another merge in the same form as the previous one
    ensemblmerge <-   GRanges(seqnames= Rle(rep(as.character(Chromosome), each=dim(snps_without_tss)[1])),
                              ranges = IRanges(start=snps_without_tss$end, end=snps_without_tss$end, names=snps_without_tss$id),
                              geneid=snps_without_tss$id, TSS=snps_without_tss$end, strand= snps_without_tss$strand)
    ASE_varsmerge <-  GRanges(seqnames=
                                Rle(rep(as.character(Chromosome), each=dim(snps_without_tss)[1])),
                              ranges = IRanges(start=snps_without_tss$end, end=snps_without_tss$end), names=snps_without_tss$id)
    #Combine the snps without tss with those of a tss to test them all later
    merge2 <- mergeByOverlaps(ASE_varsmerge,ensemblmerge)
    completemerge <- rbind (merge,merge2)
    colnames(completemerge)[2] <- "id"
    # Add the TSS and geneid information to the ASE_vars df
    newdf <- data.frame(id=completemerge$id, TSS=completemerge$TSS, Gene=completemerge$geneid)
    # This merge takes a while
    print("Now merging")
    ASE_vars <- merge(newdf, ASE_vars, by = "id")
    ASE_vars<-ASE_vars[!duplicated(ASE_vars[,c("Ind", "id")]),]
    cat("Merge complete\n")
    output_file = paste(c(output_path, "Run.model.input_Chr",  Chromosome, ".RData"), collapse="")
    print(paste(c("Saving ", output_file), collapse=""))
    
  save(list = ls(all.names = TRUE), file = output_file, envir = environment())
  print("Finished")
}  
  

