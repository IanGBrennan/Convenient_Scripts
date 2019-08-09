# in your terminal, delete the last line of each file, if there's an extra space/gap there with:
# $ sed -i '$d' PATH_TO_FILE
# or $ sed -i '$d' *.fasta # to do all fasta files

#### We'll make a loop to change the file format from Phylip to Fasta first
path = "/Users/Ian/Desktop/GenomeStripper/Elapids/Elapidae_Alignments"
out.file<-""
file.names <- dir(path, pattern =".phylip")
#dir.create(paste0(path,"/fasta")) # make a directory for the new files

extract.targets <- function(taxon, filetype){
  taxon.og <- taxon
  for (p in 1:length(file.names)) {
    # first remove the first line with the number of taxa and characters
    file <- paste(path, "/", file.names[p], sep="")
    short.file <- file.names[p]
    shortie <- strsplit(short.file, filetype)[[1]]
    #removed <- paste("sed -i '1d'", file)
    #system(removed)
    
    # read in the alignment 
    if(filetype==".fasta"){full_alignment <- read.dna(file, format="fasta", as.matrix=TRUE)}
    else if(filetype==".phylip"){taxon <- paste0(taxon.og,"\t"); full_alignment <- read.dna(file, format="sequential")}
    
    if(!taxon %in% row.names(full_alignment)){
      missing <- paste(taxon, "is not in alignment", shortie)
      print(missing)
      write.table(missing, file=paste0(path,"/",taxon,"_MISSING_LOCI.txt"), append=T, row.names=F, col.names=F)
    }
    else if(taxon %in% row.names(full_alignment)){
      # pull out just the target taxon
      target.alignment <- full_alignment[taxon,]
      
      # rename the alignment with the locus name
      row.names(target.alignment) <- shortie
      
      # write the file (appending each new sequence)
      if(filetype==".fasta") {write.FASTA(target.alignment, paste0(path, "/", taxon, "_All_Loci.fasta"), append=T)}
      else if(filetype==".phylip") {ips::write.fas(target.alignment, paste0(path, "/", taxon.og, "_All_Loci.fasta"), append=T)}
    }
  }
}

extract.targets(taxon="I12234_CTMZ_04069_Squamata_Elapidae_Pseudechis_australis__seq1", filetype=".phylip")
 