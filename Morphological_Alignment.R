library(dispRity)

## Generating a random tree
random_tree <- rcoal(10)

## Generating a random matrix
random_matrix <- sim.morpho(random_tree, characters = 50, model = "ER",
                            rates = c(rgamma, 1, 1))
## Checking the matrix scores
check.morpho(random_matrix, orig.tree = random_tree)


library(Claddis)

c.mat <- MakeMorphMatrix(random_matrix)
c.dist <- MorphDistMatrix(c.mat)
MorphMatrix2PCoA(c.dist)


library(TreeSearch); library(phangorn)
# test.mat <- ReadAsPhyDat("/Users/Ian/Desktop/Species_Tree/Reduced_Sampling/T474_Ivanov.Morphology copy.nex")
# could use the original and 'attributes(test.mat)$allLevels' to get the 'key'

# Remake the MatrixToPhyDat to inlclude a key to the labels
MatrixToPhyDat <- function (tokens) {
  allTokens <- unique(as.character(tokens))
  tokenNumbers <- seq_along(allTokens)
  names(tokenNumbers) <- allTokens
  matches <- gregexpr("[\\d\\-\\w]", allTokens, perl = TRUE)
  whichTokens <- regmatches(allTokens, matches)
  levels <- sort(unique(unlist(whichTokens)))
  whichTokens[allTokens == "?"] <- list(levels)
  contrast <- 1 * t(vapply(whichTokens, function(x) levels %in% 
                             x, logical(length(levels))))
  rownames(contrast) <- allTokens
  colnames(contrast) <- levels
  dat <- phangorn::phyDat(tokens, type = "USER", contrast = contrast)
  dat$key <- tokenNumbers
  dat
}

testo <- ReadCharacters("/Users/Ian/Desktop/GenomeStripper/Komodo_assembly/Existing_Alignments/combined_alignments/Trimmed_Alignments/Varanidae_Conrad2011_AllChars.nex")
#testo <- ReadCharacters("/Users/Ian/Desktop/GenomeStripper/Komodo_assembly/Existing_Alignments/combined_alignments/Trimmed_Alignments/Ivanov_Raw.nex")
#testo <- ReadCharacters("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/Alignments/Macropodinae_MorphOnly/Macro_Morph_Only_ALL.nex")
teste <- MatrixToPhyDat(testo)
#write.phyDat(teste, file="/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/Alignments/Macropodinae_MorphOnly/Macro_Morph_Var.nex", format="nexus", nbcol=-2)
write.phyDat(teste, file="/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/Alignments/Macropodinae_MorphOnly/Macro_Morph_Var.nex", format="nexus", nbcol=-2)


# If you want to remove some taxa before summarizing the data
remove.tips <- c("")

# Summarize the content and missing data per taxon
summarize.taxa <- function(matrix.frame = NULL, missing=c("N", "?"), taxa.threshold = 0.75){
  #lapply(matrix.frame, function(x) print(teste[x]))
  summary.list <- NULL
  # test.table <- lapply(matrix.frame, table)
  # test.toble <- lapply(test.table, function(x) names(teste$key)[match(names(test.table[[x]]), teste$key)])
  for (i in 1:length(matrix.frame)){
    char.table <- table(matrix.frame[i])
    names(char.table) <- names(matrix.frame$key)[match(names(char.table), matrix.frame$key)]
    summary.list[[i]] <- char.table
  }
  names(summary.list) <- names(matrix.frame)
  
  percentN <- NULL
  for(i in 1:length(summary.list)){
   if(missing=="N"){
     if("N" %in% names(summary.list[[i]])) {percentN[[i]] <- summary.list[[i]][names(summary.list[[i]]) == "N"] / sum(summary.list[[i]])}
     else if(! "N" %in% names(summary.list[[i]])) {percentN[[i]] <- 0}}
  else if(missing=="?") {
      if("?" %in% names(summary.list[[i]])) {percentN[[i]] <- summary.list[[i]][names(summary.list[[i]]) == "?"] / sum(summary.list[[i]])}
      else if(! "?" %in% names(summary.list[[i]])) {percentN[[i]] <- 0}}
  }
  names(percentN) <- names(matrix.frame)
  
  drop.taxa <- names(percentN[which(percentN > taxa.threshold)])

  return(list(taxon.summary = summary.list, proportion.missing = percentN, taxa2drop = drop.taxa))
}
# Summarize the content and missing data per character
summarize.characters <- function(matrix.frame = NULL, missing=c("N", "?"), character.threshold = 0.5){
  
  char.list <- NULL
  for (i in 1:(length(matrix.frame[[1]])-1)){
    char.i <- unlist(lapply(matrix.frame, '[', i))[1:(length(matrix.frame)-1)]
    char.table <- table(char.i)
    names(char.table) <- names(matrix.frame$key)[match(names(char.table), matrix.frame$key)]
    char.list[[i]] <- char.table
  }
  percentNchar <- list()
  for(i in 1:length(char.list)){
    if(missing=="N") {percentNchar[[i]] <- (char.list[[i]][names(char.list[[i]]) == "N"] / sum(char.list[[i]]))}
    else if(missing=="?") {percentNchar[[i]] <- (char.list[[i]][names(char.list[[i]]) == "?"] / sum(char.list[[i]]))}
  }
  invariables <- NULL
  for(i in 1:length(char.list)){
    inv.int <- char.list[[i]][which(!names(char.list[[i]])==missing)]
    if(length(inv.int>=2)) {invariables[[i]] <- "variable"}
    else if (length(inv.int<=1)) {invariables[[i]] <- "invariable"}
  }
  #char.list[[102]][which(!names(char.list[[102]])=="?")]
  
  drop.character <- which(percentNchar > character.threshold)
  
  return(list(proportion.missing = percentNchar, characters2drop = drop.character))
}

t.sum <- summarize.taxa(teste, missing="?", taxa.threshold = 0.75)
c.sum <- summarize.characters(teste, missing="?", character.threshold = 0.7)

# Create a new matrix removing the ill-fitting characters
remove.characters <- function(raw.object, summarized.object){
  rem.obj <- lapply(raw.object[1:length(raw.object)-1], '[', setdiff(1:length(raw.object[[1]]), summarized.object$characters2drop))
  rem.obj$key <- raw.object$key
  attributes(rem.obj)$weight <- attributes(raw.object)$weight[setdiff(1:length(raw.object[[1]]), summarized.object$characters2drop)]
  attributes(rem.obj)$class <- "phyDat"; attributes(rem.obj)$nr <- length(rem.obj[[1]]); 
  attributes(rem.obj)$levels <- attributes(raw.object)$levels; attributes(rem.obj)$allLevels <- attributes(raw.object)$allLevels
  attributes(rem.obj)$type <- "User"; attributes(rem.obj)$contrast <- attributes(raw.object)$contrast
  return(rem.obj)
}

teste3 <- remove.characters(teste, c.sum)

summarize.taxa(teste3, missing="?", taxa.threshold = 0.7)
      teste3 <- teste3[!names(teste3)=="key"]

write.phyDat(teste3, file="/Users/Ian/Desktop/GenomeStripper/Komodo_assembly/Existing_Alignments/combined_alignments/Trimmed_Alignments/Varanidae_Conrad_Reduced.nex", format="nexus", nbcol=-2)

sum.test <- summarize.matrix(teste)


char1 <- unlist(lapply(teste, '[[', 2))[1:(length(teste)-1)]
names(char1) <- NULL

unlist(teste[[i]])[1]


sum(sum.test[[2]])

######################################################
## Come back to this because it's probably valuable
######################################################
testi <- ReadMorphNexus("/Users/Ian/Desktop/GenomeStripper/Komodo_assembly/Existing_Alignments/combined_alignments/Trimmed_Alignments/Ivanov_Raw.nex")

testid <- lapply(test.mat, "2" <- "?")

test.phydat <- as.phyDat(testo, type="USER", levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "-"), ambiguity = c("?", "N"))

meb <- test.mat$Shinisaurus_crocodilurus
meb[which(meb==2)] <- "?"
