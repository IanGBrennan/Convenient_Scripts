# function will split up the tip labels from AHE trees into usable info
# save it externally as a CSV for later!

extract.labels <- function(phy) {
  labels <- phy$tip.label
  data.table <- NULL
  for (i in 1:Ntip(phy)) {
    int.mat <- NULL
    splitup <- strsplit(labels[i], split="_")
    ahe.index <- splitup[[1]][1]; int.mat[1] <- ahe.index
    museum.index <- splitup[[1]][2]; int.mat[2] <- museum.index
    order <- splitup[[1]][3]; int.mat[3] <- order
    family <- splitup[[1]][4]; int.mat[4] <- family
    genus <- splitup[[1]][5]; int.mat[5] <- genus
    species <- splitup[[1]][6]; int.mat[6] <- species
    genus_species<- paste(genus, species, sep="_"); int.mat[7] <- genus_species
    int.mat[8] <- labels[i]
    data.table <- rbind(data.table, t(as.data.frame(int.mat)))
  }
  colnames(data.table) <- c("ahe_index", "museum_index", "order", "family", 
                            "genus", "species", "original_genus_species", "tip_label")
  rownames(data.table) <- NULL
  return(data.table)
}

extract.labels.SqCL <- function(phy) {
  labels <- phy$tip.label
  data.table <- NULL
  for (i in 1:Ntip(phy)) {
    int.mat <- NULL
    splitup <- strsplit(labels[i], split="_")
    #ahe.index <- splitup[[1]][1]; int.mat[1] <- ahe.index
    #museum.index <- splitup[[1]][2]; int.mat[2] <- museum.index
    #order <- splitup[[1]][3]; int.mat[3] <- order
    family <- splitup[[1]][1]; int.mat[1] <- family
    genus <- splitup[[1]][2]; int.mat[2] <- genus
    species <- splitup[[1]][3]; int.mat[3] <- species
    genus_species<- paste(genus, species, sep="_"); int.mat[4] <- genus_species
    museum.id <- splitup[[1]][4]; int.mat[5] <- museum.id
    museum.no <- splitup[[1]][5]; int.mat[6] <- museum.no
    museum.index <- paste(museum.id, museum.no, sep="_"); int.mat[7] <- museum.index
    int.mat[8] <- labels[i]
    data.table <- rbind(data.table, t(as.data.frame(int.mat)))
  }
  colnames(data.table) <- c("Family", "Genus", "Species", "Genus_species", "Museum_ID", 
                            "Museum_Number", "Museum_Index", "Name_in_Tree")
  rownames(data.table) <- NULL
  return(data.table)
}

extract.labels.basic <- function(phy) {
  labels <- phy$tip.label
  data.table <- NULL
  for (i in 1:Ntip(phy)) {
    int.mat <- NULL
    splitup <- strsplit(labels[i], split="_")
    #ahe.index <- splitup[[1]][1]; int.mat[1] <- ahe.index
    #museum.index <- splitup[[1]][2]; int.mat[2] <- museum.index
    #order <- splitup[[1]][3]; int.mat[3] <- order
    #family <- splitup[[1]][4]; int.mat[4] <- family
    genus <- splitup[[1]][1]; int.mat[1] <- genus
    species <- splitup[[1]][2]; int.mat[2] <- species
    genus_species<- paste(genus, species, sep="_"); int.mat[3] <- genus_species
    museum.id <- "NA"; int.mat[4] <- museum.id
    museum.no <- splitup[[1]][3]; int.mat[5] <- museum.no
    museum.index <- paste(museum.id, museum.no, sep="_"); int.mat[6] <- museum.index
    int.mat[7] <- labels[i]
    data.table <- rbind(data.table, t(as.data.frame(int.mat)))
  }
  colnames(data.table) <- c("Genus", "Species", "Genus_species", "Museum_ID", 
                            "Museum_Number", "Museum_Index", "Name_in_Tree")
  rownames(data.table) <- NULL
  return(data.table)
}

# write.csv('output', file='...')