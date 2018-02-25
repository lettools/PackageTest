# Package Test

`PackageTest` is a package containing two functions that calculates statistical associations between the level of expression at ASE sites and the presence of *cis* sequence variants by applying running permutations of the input data and applying a generalized linear model (GLM) each time.


### Dependencies ###
* This package was built in R/3.2.2 and requires the following packages to run:
* library(GenomicRanges)
* library(reshape2)
* library(rtracklayer)
* library(speedglm)
* library(stringr)
* library(plyr)
* library(biomaRt)

### Input files 
Only those columns that are required are shown. The sample, legend and hap files are in the IMPUTE format (https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html).

Both the legend and sample files have headers, but the hap does not. Your hap file should have the same number of rows as that of the legend (once the legend's header is accounted for) and should also have twice the number of columns as the sample has rows, because each individual will have two haplotypes. Below are examples of acceptable formats. 

The ASE file must contain the position of the ASE sites, including the chromosome they are on, the individuals in which they are found, and the alleles at the different haplotypes, including their respective read counts.

```
# Format of sample Samples file
#  ID        SEX    POP
# HG00096    male    GBR
# HG00097  female    GBR
# HG00099  female    GBR
# NA20827    male    TSI


# Format of the legend file
#                   id     position    a0  a1   TYPE
#           10:60515:C:T    60515       C   T   Biallelic_SNP
#  rs148087467:60523:T:G    60523       T   G   Biallelic_SNP
# rs147855157:61372:CA:C    61372       CA  C   Biallelic_INDEL


# Format of the hap file
#0 0 1 0 0 0 0 0 0 0
#0 0 0 0 1 0 0 1 0 0
#0 0 0 0 0 0 0 0 0 0


# Format of sample ASE file
#   chr    end    ref    alt      Ind         refCount    altCount
#  22   135032   G      A       HG00276       19            0
#  22   135032   G      A       HG00282       12            0
#  22   135032   G      A       NA11831       10            0
```


### Functions
1. `Gen.input`is a function to generate an input for the model to run. Separating out into two functions saves memory and greatly speeds up parallelisation for what is a computationally demanding task in the second function. The input files are formatted, gene and transcript start site information added to the ASE file, and finally output in .RData format to be loaded into the second function.
2. `Run.Model` is a function that allows you to measure statistical association between nearby regulatory variants and the level of expression at a heterozygous coding polymorphism, controlling for factors such as sex and population, by utilising a generalized linear model and applying permutations to the data in order to provide a robust p-value. As the function supports parallelisation, a number of .txt files equal to the number of tasks, `numTasks`, specified in the function will be outputted.  


## Installing in R ##
```r
# First, load the devtools library: 
library(devtools)

# Then, install the package from GitHub:
install_github("ac1990/PackageTest")

# Then load package as normal:
library("PackageTest")

# Load the sample files 
download.file("https://www.dropbox.com/s/cjb7xsx0qvos0y4/sample_ASE.txt.gz?dl=1","sample_ASE.txt.gz")
download.file("https://www.dropbox.com/s/dp173n4af6fo6c8/sample_haplotypes.hap.gz?dl=1","sample_haplotypes.hap.gz")
download.file("https://www.dropbox.com/s/9ngy0ao6hahpic2/sample_Samples.txt.gz?dl=1","sample_Samples.txt.gz")
download.file("https://www.dropbox.com/s/w2arx4728c3rycv/sample_legend.leg.gz?dl=1","sample_legend.leg.gz")


# Run the first function, Gen.input, with the sample files 
# This function allows you to output .RData files to be used as input for the analysis run in the second function.
# This only needs to be done once for each chromosome. This ensures the running of the model is not too computationally demanding.
# Parameters such as the number of permutations and the TSS window can be tweaked in the Run.model function. 
# The processing of input files need not be done repeatedly, as they are constant
# Hence separating this processing out allows tweaks to be made to the parameters without having to re-run the processing.
# Usage:
     Gen.input(Chromosome, ASE_file, legend_file, haplotypes_file, Samples,
       output_path, Species = "hsapiens", EnsemblVersion = NULL)
       
# Arguments:
# Chromosome: Specify chromosome. Must be in the same format as that found in Ensembl. For example, if you are
# looking at chromosome 22, it would either be `Chromosome = 22` or simply `22`. If looking at the x chromosome,
# it would be `Chromosome="x"` or `"x"`
# ASE_file: Specify file containing allele-specific expression information
# legend_file: Specify legend file containing SNP IDs and position info etc
# haplotypes_file: Provide haplotypes corresponding to SNPs in legend file, for given individuals found in samples file
# Samples: File containing cohort information
# output_path: Directory path to save output files to
# Species: Select same species as cohort. Defaults to "hsapiens"
# EnsemblVersion: Specify version of Ensembl to download from. Defaults to NULL which is the most recent build

# Examples:

# Downloading Ensembl information for chromosome 22 of hsapiens, using the build corresponding to the sample data
     Gen.input(22, "path_to_ASE_file/ASEfile.txt.gz",
                         "path_to_legend_file/file.legendfile.gz",
                        "path_to_haplotypes_file/hapfile.hap.gz",
                         "path_to_samples_file/samples.txt", EnsemblVersion=78)
     
# Downloading Ensembl information for chromosome 22 of hsapiens, using the most recent build 
     Gen.input(22, "path_to_ASE_file/ASEfile.txt.gz",
                         "path_to_legend_file/file.legendfile.gz",
                        "path_to_haplotypes_file/hapfile.hap.gz",
                         "path_to_samples_file/samples.txt")
     
# Downloading Ensembl information for chromosome 1 of mmusculus, using the most recent build
     Gen.input(1, "path_to_ASE_file/ASEfile.txt.gz",
                         "path_to_legend_file/file.legendfile.gz",
                       "path_to_haplotypes_file/hapfile.hap.gz",
                         "path_to_samples_file/samples.txt",



# Next, run the model.
# The RData files outputted from the function, Gen.input, can now be used to run the analysis in the second function, below.
# Usage:
 Run.Model(input_file, Task, progress_path, output_path, numTasks = 1000,
       Chromosome, numPerms = 100000, TSSwindow = 500000, pval_threshold = 0.00005)   
             
# Arguments: 
# input_file: This is the RData file outputted from the first function, Gen.input

#    Task: This analysis is very computationally burdensome. To speed up the process it is an advantage to split 
#          it up into tasks that may be run on multiple nodes, concurrently.

# progress_path: The function allows checkpointing, which allows the
#                program to be killed and picked up again at a later stage,
#                without starting from the beginning. Checkpointing occurs
#                every 2 hours

# output_path: Specify path to save the results to

# numTasks: Select how many jobs/tasks to split the process up into.
#           Defaults to 1000

# Chromosome: Specify chromosome used in input file

# numPerms: Select how many permutations to run. Along with splitting the
#          process up into simultaneous tasks, this is the biggest
#         factor in determining how long the analysis will take.
#          However, the more permutations, in general and up until a
#          point, the more precise and accurate the results may be; for
#          example, if set to 100, the minimum p-value that can possibly
#          be reached as a result of permutations, is 0.01. Defaults to
#          100,000.

# TSSwindow: This represents the distance over which nearby variants will
#          be selected, either side of the transcript start site.
#          Defaults to 500kb

# pval_threshold: There is a theoretical minimum p-value for each
#                 particular combination of reference and alternative alleles
#                 for a given set of individuals for a given nearby variant of
#                 an ASE site.  This parameter sets the upper limit. Default is
#                 0.00005. In this example the model will not be run if it is
#                 not possible to reach a p-value as low as 0.00005, even
#                 theoretically.
          
          
# Examples:          
          
# Run model with task set to 10, chromosome to 22, for 100,000 permutations, a transcript start site window of 500kb and a theoretical p-value threshold of 0.00005
     Run.Model("input_file.RData", 10, 
             "path_to_progress_file",
             "output_path",
             Chromosome=22)
     
# Run model with task    set to 10, chromosome to 22, for 10,000 permutations, a transcript start site window of 500kb and a theoretical p-value of 0.0005
     Run.Model("input_file.RData", 10,
             "path_to_progress_file",
             "output_path",
             Chromosome=22, numPerms=10000, pval_threshold=10000)
     
# Run model with task    set to 2, chromosome to 12, for 10,000 permutations, a transcript start site window of 1Mb and a theoretical p-value of 0.0005
     Run.Model("input_file.RData", 2,
             "path_to_progress_file",
             "output_path",
             Chromosome=22, numPerms=10000,TSSwindow=100000, pval_threshold=10000)
     
#     NB: The smallest possible p-value attainable as a result of running permutations is 1/numPerms. 
#         Hence, there is no advantage to setting the minimum p-value threshold to below this number.
```



## Output format ##

The first four columns are the ID of the ASE site, its position, the individual, and the ID of the nearby variant. Then follows 7 columns representing the test statistics of the intercept, reads, sex, sub-populations and variant, with the next 7 columns representing their corresponding p-values. The next two columns are the transcript start site of the ASE site and the gene to which the nearby variant belongs. Where the nearby variant is not found in a gene, the ID of the variant itself is given. The final two columns are the number of permutations, and the number of permutations in which the test statistic is exceeded.  

### Selecting significant sites ###
To calculate the permuted p-values, in R do the following to your results dataframe, `df`:
```
# Divide numPermExceed column by the numPerm column
df2 <- df[which(as.numeric(df$numPerm) > 0),]
df2$permuted_p <-as.numeric(df2$numPermExceed)/as.numeric(df2$numPerm)
```
results <- results













## Plotting function to be done later ##
```
# Stick this in plotting function maybe
# To obtain false discovery rate (fdr) adjusted p-values:
df2$adjusted_p <- p.adjust(df2$permuted_p,"fdr",n=length(df2$permuted_p))
```




##### Things to set as options:

9. If you want to look at the nearby allele on the same chr as coding het site or other chr
10. Also, if they don't provide info like POP and SEX they can still do the analysis, so need to modify the glm equation
11. The ability to turn off checkpointing



Could even put these all in a function. They call the function and it installs them all.




EOF
