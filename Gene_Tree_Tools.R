library(ape); library(phytools)

files <- list.files(pattern=".trees", recursive=FALSE)

setwd("/Users/Ian/Google.Drive/ANU Herp Work/Pygopodidae_ST/Pygopodidae.Species.Files/starBEAST2/Nuclear_Chain3")

# remember to set your working directory before using the function!
batch_consensus <- function(file.extension=".trees", burnin.percent, 
                            path.to.treeannotator="/Applications/BEAST_2.5.0/bin/treeannotator") {
  files <- list.files(pattern=file.extension, recursive=FALSE)
  
  for(i in 1:length(files)){
    input = files[i]
    #in BEAST2, this is a percentage, might need to adjust this script if using BEAST1 
    exe = path.to.treeannotator
    #this is the location of treeannotator in your current BEAST directory
    filename <- strsplit(input, file.extension)[[1]]
    output <- paste(filename, "con.tre", sep=".")
    burn.in <- paste("-b",burnin.percent)
    call <- paste(exe, burn.in, input, output) 
    system(call)
  }
}

batch_consensus(file.extension=".trees", burnin.percent=10)
