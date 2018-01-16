

#### We'll make a loop to randomly sample 10 loci 100 times
setwd("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Protea_Files/OnePer_Sequence_Files/Actual_alignments")
path = "/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Protea_Files/OnePer_Sequence_Files/Actual_alignments"
out.file<-""
file.names <- dir(path, pattern=".phylip")
for (q in 1:100) {
  pick10 <- paste(sample(file.names, 10, replace=F), collapse=" ")
  #pick <- paste(shQuote(pick10), collapse=",")
  outname <- paste("Protea_random10_num_",q,".phy", sep="")
  call <- paste("cat", pick10, ">>", outname)
  system(call)
}

tree <- read.tree("/Applications/MCMCtree_tutorial/mcmctree_Protea/Protea_astral_input.tre")
tree[[102]]$edge.length <- NULL
write.tree(tree[[102]], "/Applications/MCMCtree_tutorial/mcmctree_Protea/Protea_input.tre")
test <- tree[[102]]
sort(test$tip.label)
test <- drop.tip(test, c("_Protea_nitida_seq1_1447", "rotea_rupicola_seq1_161"))
write