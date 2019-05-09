# We want to clean up a bunch of filthy alignments
# They have taxa with no sequence data, lots of 
# Columns with missing data, and etc.

library(ape); library(ips)

# set the working directory
setwd("/Users/Ian/Desktop/Gecko_AHE/Original_Single_line_FASTA/Cleaned_Alignments")

# Create a function like "deleteGaps" to drop columns of missing data (- or N) with a threshold
deleteWeirdos <- function(x, nmax=nrow(x)-4){
  if (!inherits(x, "DNAbin")) 
    stop("'x' is not of class 'DNAbin'")
  gapCounter <- function(n, nmax) {
    length(which(n %in% as.raw(240) | n %in% as.raw(4))) > nmax
  }
  id <- apply(x, 2, gapCounter, nmax = nmax)
  x[, !id]
}
# this function works because of these characters:
# as.integer(as.DNAbin("-")) # this is the binary coding for a gap "-"
# as.integer(as.DNAbin("N")) # this is the binary coding for missing base "N"


file.extension <- ".fasta"
files <- dir(getwd(), pattern=file.extension)

for (p in 1:length(files)){
  
  # get the locus name
  locus.name <- str_split(files[p], pattern = file.extension)[[1]][1]
  
  # read in the alignment file
  aligned <- read.dna(files[p], format="fasta", as.matrix=TRUE)
  
  # remove empty rows and columns
  nonempty <- deleteEmptyCells(aligned, quiet=T)
  
  # make sure there's still sequences:
  if (length(nonempty) > 0) {
    
    # clean up columns that are primarily missing data (- and N)
    cleaned <- deleteWeirdos(nonempty, nmax = round(nrow(nonempty)/2))
    
    # write the file (appending each new sequence)
    write.FASTA(cleaned, paste0("Clean_", files[p]))
    
    # I started desiging a partition file, but easier to do with 'concat' in AMAS
    #part.info <- paste("charset", locus.name, "=", )
    
  } else (print(paste("alignment", files[p], "had no sequence data")))
  
  print(paste(files[p],";", (p+1)))
}




###################################################################
# You could do this less efficiently, but more flexibly this way,
# it should give you a similar resulting alignment
files <- dir(getwd(), pattern=".fasta")

for (p in 1:length(files)){
  
  # read in the alignment file
  aligned <- read.dna(files[p], format="fasta", as.matrix=TRUE)
  
  # trim off hanging ends
  trimmed <- trimEnds(aligned, min.n.seq = round(nrow(aligned)/2))
  
  # remove gappy bits
  gapped <- deleteGaps(trimmed, nmax=round(nrow(tri)/2))
  
  # write the file (appending each new sequence)
  write.FASTA(gapped, paste0("TrimmedGapped_", files[p]))
  
  # view the alignment
  # image(gapped)
  
  print(paste(files[p],";", (p+1)))
}


## Trim out reference sequences
align.path <- "/Users/Ian/Desktop/Gecko_AHE/Clean_FASTAs/Lialis_sequences"
file.names <- dir(align.path, pattern =".fasta")
output.dir <- "Lialis_burtonis_sequences"
dir.create(paste0(align.path,"/", output.dir)) # make a directory for the new files

remove.ref <- function(taxon, filetype){
  for (p in 1:length(file.names)) {
    # first remove the first line with the number of taxa and characters
    file <- paste(align.path, "/", file.names[p], sep="")
    short.file <- file.names[p]
    shortie <- strsplit(short.file, filetype)[[1]]
    #removed <- paste("sed -i '1d'", file)
    #system(removed)
    
    # read in the alignment 
    full_alignment <- read.dna(file, format="fasta", as.matrix=TRUE)
    
    
    if(!taxon %in% row.names(full_alignment)){
      missing <- paste(taxon, "is not in alignment", shortie)
      print(missing)
      write.table(missing, file=paste0(align.path,"/",taxon,"_MISSING_LOCI.txt"), append=T, row.names=F, col.names=F)
    } else if(taxon %in% row.names(full_alignment)){
      # pull out just the target taxon
      target.alignment <- full_alignment[taxon,]
      
      # rename the alignment with the locus name
      row.names(target.alignment) <- shortie
      
      # write the file (appending each new sequence)
      write.FASTA(target.alignment, paste0(align.path,"/",output.dir, "/", taxon, "_", shortie, ".fasta"))
      print(paste("isolated alignment", p, "!"))
    }
  }
}

remove.ref(taxon="Lialis_burtonis", filetype=".fasta")


