source("https://bioconductor.org/biocLite.R")
biocLite("muscle")
library(muscle)
biocLite("Biostrings")
library(Biostrings)
library(ape)
library(BioGeoBEARS)
library(phylotools)

#### This preliminary step is for eliminating 'duplicate' sequences from phased data
# doing this to trim the python data down to a single sequence per individual, to agree with the elapids
path = "/Users/Ian/Desktop/combined/to.drop" # make sure to designate you path appropriately!
file.names <- dir(path, pattern =".phylip")
for (m in 1:length(file.names)) {
  file <- paste(path, "/", file.names[m], sep="")
  drop.one <- paste("sed -i '' 'n; d'", file) # this drops even lines, drop odd lines with "sed '1d; n; d'"
  system(drop.one)
}

#### Start by combining the loci that are orthologous, into a single Phylip file
overlap <- read.csv("loci.overlap.csv", header=T) # read in the congruence file
for (z in 1:nrow(overlap)) {
  locus.a <- overlap[z,1]
  locus.b <- overlap[z,2]
  
  file.a <- paste(locus.a, "_phylip.txt", sep="")
  file.b <- paste(locus.b, ".phylip", sep="")
  
  new.name <- paste("combined.",locus.a,"__",locus.b,".fasta", sep="")
  call <- paste("cat", file.a, file.b, ">", new.name)
  system(call)
}

# ignore this, just a method for reading in Phylip into R
zed <- read.phylip("T222_L100_phylip.txt")
head(zed)

#### We'll make a loop to change the file format from Phylip to Fasta first
path = "/Users/Ian/Desktop/combined"
out.file<-""
file.names <- dir(path, pattern =".phylip")
for (q in 4:length(file.names)) {
  # first remove the first line with the number of taxa and characters
  file <- paste(path, "/", file.names[q], sep="")
  short.file <- file.names[q]
  removed <- paste("sed -i '' '1d'", file)
  system(removed)
  
  #alternatively this can be done more gracefully with 'system2'
  #system2("sed", args=c("-i", "''", "'1d'", file))
  
  #now we'll add the 'fasta carrot'
  arrow <- paste("sed -i '' 's/^/>/'", file)
  system(arrow)
  
  #then insert a return for each tab and save the file with '.fasta' ending
  fasta.name <- paste(short.file, ".fasta", sep="")
  returned <- paste("tr '\t' '\n' <", file, ">", fasta.name)
  system(returned)
}

#### Now we'll loop through sending each fasta file to Muscle to be aligned (this could take a while)
fasta.files <- dir(path, pattern=".fasta")
for(i in 1:length(fasta.files)) {
  locus <- paste(path,"/",fasta.files[1], sep="")
  new.name <- paste("aligned.", fasta.files[1], sep="")
  call <- paste("/Users/Ian/Desktop/combined/muscle3.8.31_i86darwin64",
                "-in", locus, "-out", new.name)
  system(call)
}

