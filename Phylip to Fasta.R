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