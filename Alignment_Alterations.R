name.changes <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Cichlids/Cichlid_name_changes.csv",
                         header = T)

setwd("/Users/Ian/Desktop/Renamed_415_alignments") # set your working directory

# read in all the files you want to handle
#(designate the folder, then isolate the files by their standard ending)
files <- list.files(pattern=".phy", recursive=FALSE)

#################################################################
# loop through all the alignment files in the folder
# use sed to remove the prefixes on the names
#################################################################
path = "/Users/Ian/Desktop/Renamed_415_alignments"
out.file <- ""
file.names <- dir(path, pattern = ".phy")
for(i in 1:length(file.names)){
  in.file <- paste(path, "/", file.names[i], sep="")
  input = file.names[i]

  filename <- tail(unlist(strsplit(input, ".phy")), 1)
  datapath <- gsub(filename, "", input)
  filecall <- head(unlist(strsplit(input, "ENSONIE0000")), 2)[2]
  #output <- paste(path, "/", "Locus_", filecall, sep="")

  pre.call <- paste("sed -i '' 's/^[^'",filename, "_","']*'",filename, "_","'//'", sep="")

  call <- paste(pre.call, in.file)
  system(call)
}

#################################################################
# loop through all the alignment files in the folder
# use unix to rename all the file names
#################################################################
for(y in 1:length(file.names)) {
  #in.file <- paste(path, "/", file.names[y], sep="")
  input = file.names[y]
  
  filename <- tail(unlist(strsplit(input, ".phy")), 1)
  datapath <- gsub(filename, "", input)
  filecall <- head(unlist(strsplit(input, "ENSONIE0000")), 2)[2]
  #output <- paste(path, "/", "Locus_", filecall, sep="")
  output <- paste("Locus_", filecall, sep="")
  
  rename <- paste("mv", input, output)
  system(rename)
  
  #output <- paste("Locus_", filecall, sep="")
  #renamecall <- paste("rename -f 's/", input, "/", output,"'")
  #system(renamecall)
}

#################################################################
# loop through all the alignment files in the folder
# use perl to change the taxon names according to a CSV file (name.changes)
#################################################################
setwd("/Users/Ian/Desktop/Goodnames_415_alignments") # set your working directory
for(k in 1:nrow(name.changes)){
  current <- name.changes[k,]
  old.name <- current[,1]
  new.name <- current[,2]
  
  rename <- paste("perl -pi -w -e 's/", old.name, "/", new.name, "/g;' *.phy", sep="")
  system(rename)
}

#################################################################
# loop through all the alignment files in the folder
# use perl to drop some taxa according to a CSV file (name.changes)
#################################################################
setwd("/Users/Ian/Desktop/Reduced_415_alignments") # set your working directory
to.drop <- name.changes[,3][1:78]
for(j in 1:length(to.drop)){
  drop.it <- to.drop[j]

  drop.sample <- paste("sed -i '' -e '/", drop.it, "/d' *.phy", sep="")
  system(drop.sample)
}
### when you're done, remember to change the number of samples in each alignment
#### at the top of your phylip files!
setwd("/Users/Ian/Desktop/Reduced_415_alignments") # set your working directory
call <- paste("perl -pi -w -e 's/ 147 / 69 /g;' *.phy")
system(call)

#################################################################
# finally, create a RAxML gene tree for each locus
# by calling RAxML for each alignment (change settings as you see fit)
#################################################################
setwd("/Users/Ian/Desktop/Reduced_415_alignments") # set your working directory
alignments <- list.files(pattern=".phy", recursive=FALSE)

for (z in 84:length(alignments)) {
  cat("tree", z, "of", length(alignments), "\n") #keep track of what tree/loop# we're on
  current <- alignments[z]
  locus.name <- tail(unlist(strsplit(current, ".phy")), 1)
  
  tree <- paste("/Users/Ian/Desktop/raxmlHPC-PTHREADS-SSE3  -f a -T 8 -m GTRGAMMA -s ", current, " -x 12345 -N 100 -k -n ", locus.name)
  system(tree)
}
