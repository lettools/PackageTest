#example sample data
ID <- c("HG00096", "HG00097", "HG00099", "NA20827") 
SEX <- c("male", "female", "female", "male")
POP <- c("GBR", "GBR", "GBR", "TSI")
sample_example <- data.frame(ID, SEX, POP) 

#example allele counts data
chr <- c(22,22,22)
start <- c(135031, 135031, 135031)
end <- c(135032, 135032, 135032)
ref <- c("G","G","G")
alt <- c("A","A","A")
Ind <- c("HG00276", "HG00282", "NA11831")
refCount <- c(19,12,10)
altCount <- c(0,0,0)
ASE_example <- data.frame(chr, start, end, ref, alt, Ind, refCount, altCount)

#example legend data
id <- c("10:60515:C:T", "rs148087467:60523:T:G", "rs147855157:61372:CA:C")
position <- c(60515, 60523, 61372)
a0 <- c("C", "T", "CA")
a1 <- c("T", "G", "C")
TYPE <- c("Biallelic_SNP", "Biallelic_SNP", "Biallelic_INDEL")
LEG_example <- data.frame(id, position, a0, a1, TYPE)
