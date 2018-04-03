library(sp)
library(adehabitatHR)
library(rgeos)
library(rworldmap); library(ggmap)
library(phangorn)
library(R.utils)
library(dplyr)

#data(Anolis.data)
##Create a geography.object with a modified edge matrix
##First, specify which region each branch belonged to:
#Anolis.regions<-c(rep("cuba",14),rep("hispaniola",17),"puerto_rico")
#map<-cbind(Anolis.data$phylo$edge,Anolis.regions)
#  # or  #
#map <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/ANOLIS_TEST_DATA.csv", header=T)
#    map <- map[,c("Name_in_Tree", "Latitude", "Longitude")]
##CreateGeoObject(Anolis.data$phylo,map=Anolis.map)
#phylo <- Anolis.data$phylo
#
#distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/ANOLIS_TEST_DATA.csv", header=T)
#    distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]
#        distribution <- distribution[complete.cases(distribution),] 
#
CreateGeoObject_SP <- function (phylo, map) {
  if (any(grepl("___", phylo$tip.label))) {
    stop("script will not work with '___' in tip labels; remove extra underscores")
  }
  TWOLETTERS <- paste(rep(LETTERS, each = 26), LETTERS, sep = "")
  THREELETTERS <- paste(rep(TWOLETTERS, each = 26), LETTERS, 
                        sep = "")
  nodeDist <- vector(mode = "numeric", length = phylo$Nnode) # vector of the internal nodes
  totlen <- length(phylo$tip.label) # the number of tips
  root <- totlen + 1 # number of tips +1 for the root
  heights <- nodeHeights(phylo)
  for (i in 1:dim(phylo$edge)[1]) {
    nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i] # heights above root of each internal node
  }
  nodeDist <- c(nodeDist, max(heights)) # add the total height on the end of the list of nodeheights
  nodeDiff <- diff(nodeDist) # get the distances between the internal nodes
  flag = 0
  if (sum(nodeDiff < 0) > 0) {
    node.order <- match(rank(heights[, 1], ties.method = "min"), 
                        seq(1, by = 2, len = phylo$Nnode))
    node.order <- node.order + totlen
    old.edge <- phylo$edge
    old.phylo <- phylo
    phylo$edge[, 1] <- node.order
    for (j in 1:length(phylo$edge[, 2])) {
      if (phylo$edge[j, 2] > totlen) {
        phylo$edge[j, 2] <- phylo$edge[, 1][match(phylo$edge[j, 
                                                             2], old.edge[, 1])]
      }
    }
    nodeDist <- vector()
    for (i in 1:dim(phylo$edge)[1]) {
      nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i]
    }
    nodeDist <- c(nodeDist, max(heights))
    nodeDiff <- diff(nodeDist)
    flag = 1
  }
  mat <- matrix(nrow = 0, ncol = 3)
  counter_three_letters <- 0
  
  # I can copy the names for MRCAs here and apply them to nodes on the tree 
  for (i in 1:phylo$Nnode) {
    other <- phylo$edge[phylo$edge[, 1] == i + totlen, 2]
    for (b in other) {
      int <- matrix(ncol = 3)
      int[1] <- i + totlen
      if (b > totlen) {
        counter_three_letters <- counter_three_letters + 
          1
        int[2] <- paste(".", THREELETTERS[counter_three_letters], 
                        sep = "")
        int[3] <- b
      }
      else {
        int[2] <- phylo$tip.label[[b]]
        int[3] <- 0
      }
      mat <- rbind(mat, int)
    }
  }
  #if (any(class(map) == "phylo")) {
  #  for (i in 1:length(mat[, 1])) {
  #    if (mat[i, 3] == 0) {
  #      mat[i, 3] <- as.character(match(mat[i, 2], phylo$tip.label))
  #    }
  #  }
  #  maps.object <- map$maps
  #  brchange <- vapply(maps.object, function(x) any(length(names(x)) > 
  #                                                    1), 1)
  #  if (sum(brchange) != 0) {
  #    newtimes <- vector()
  #    for (i in 1:length(brchange)) {
  #      if (brchange[i] != 0) {
  #        intlen <- length(maps.object[[i]])
  #        inttime <- vector()
  #        for (j in 1:(intlen - 1)) {
  #          inttime <- c(inttime, maps.object[[i]][j])
  #          nt <- as.numeric(heights[i, 1] + sum(inttime))
  #          newtimes <- c(newtimes, nt)
  #        }
  #      }
  #    }
  #    old.Dist <- nodeDist
  #    old.Diff <- nodeDiff
  #    nodeDist <- sort(c(nodeDist, newtimes))
  #    nodeDiff <- diff(nodeDist)
  #  }
  #  nat <- list()
  #  nodecount = 1
  #  for (i in 1:length(nodeDiff)) {
  #    if (i == 1) {
  #      hold.m <- mat[as.numeric(mat[, 1]) <= (totlen + 
  #                                               nodecount), c(2, 3)]
  #      int <- dim(hold.m)[1]
  #      for (m in 1:int) {
  #        hold.m[m, 2] <- names(maps.object[[which(phylo$edge[, 
  #                                                            2] == as.numeric(hold.m[m, 2]))]])[1]
  #      }
  #      nat[[i]] <- hold.m
  #    }
  #    else {
  #      if (nodeDist[i] %in% old.Dist) {
  #        P <- mat[as.numeric(mat[, 1]) <= (totlen + 
  #                                            nodecount), c(2, 3)]
  #        hold.m <- rbind(P[as.numeric(P[, 2]) <= totlen, 
  #                          ], P[as.numeric(P[, 2]) > (totlen + nodecount), 
  #                               ])
  #        int <- dim(hold.m)[1]
  #        for (m in 1:int) {
  #          iden <- which(phylo$edge[, 2] == as.numeric(hold.m[m, 
  #                                                             2]))
  #          if (brchange[iden] == 0 || !(hold.m[m, 1] %in% 
  #                                       nat[[i - 1]][, 1])) {
  #            hold.m[m, 2] <- names(maps.object[[iden]])[1]
  #          }
  #          else {
  #            num = 1
  #            while (round(nodeDist[i + 1] - old.Dist[phylo$edge[iden, 
  #                                                               1] - totlen], 6) > round(sum(maps.object[[iden]][1:num]), 
  #                                                                                        6)) {
  #              num = num + 1
  #            }
  #            hold.m[m, 2] <- names(maps.object[[iden]])[num]
  #          }
  #        }
  #        nat[[i]] <- hold.m
  #      }
  #      else {
  #        P <- mat[as.numeric(mat[, 1]) <= (totlen + 
  #                                            nodecount), c(2, 3)]
  #        hold.m <- rbind(P[as.numeric(P[, 2]) <= totlen, 
  #                          ], P[as.numeric(P[, 2]) > (totlen + nodecount), 
  #                               ])
  #        int <- dim(hold.m)[1]
  #        for (m in 1:int) {
  #          iden <- which(phylo$edge[, 2] == as.numeric(hold.m[m, 
  #                                                             2]))
  #          if (brchange[iden] == 0 | !(hold.m[m, 1] %in% 
  #                                      nat[[i - 1]][, 1])) {
  #            hold.m[m, 2] <- names(maps.object[[iden]])[1]
  #          }
  #          else {
  #            num = 1
  #            while (round(nodeDist[i + 1] - old.Dist[phylo$edge[iden, 
  #                                                               1] - totlen], 6) > round(sum(maps.object[[iden]][1:num]), 
  #                                                                                        6)) {
  #              num = num + 1
  #            }
  #            hold.m[m, 2] <- names(maps.object[[iden]])[num]
  #          }
  #        }
  #        nat[[i]] <- hold.m
  #      }
  #    }
  #    if (nodeDist[i + 1] %in% old.Dist) {
  #      nodecount = nodecount + 1
  #    }
  #  }
  #  geography.matrix <- list()
  #  for (i in 1:length(nodeDiff)) {
  #    timei <- unlist(nat[[i]])
  #    var.list <- timei[, 1]
  #    len = length(var.list)
  #    int.mat <- matrix(nrow = len, ncol = len)
  #    rownames(int.mat) <- var.list
  #    colnames(int.mat) <- var.list
  #    diag(int.mat) <- 1
  #    for (j in 1:length(var.list)) {
  #      for (k in 1:length(var.list)) {
  #        if ((lower.tri(int.mat, diag = TRUE)[j, k] == 
  #             TRUE)) {
  #          int.mat[j, k] <- ifelse(timei[match(var.list[j], 
  #                                              timei[, 1]), 2] == timei[match(var.list[k], 
  #                                                                             timei[, 1]), 2], 1, 0)
  #        }
  #      }
  #    }
  #    int.mat[upper.tri(int.mat) == TRUE] <- t(int.mat)[upper.tri(t(int.mat)) == 
  #                                                        TRUE]
  #    geography.matrix[[i]] <- int.mat
  #  }
  #}
  if (is.data.frame(map)) {
    nat <- list()
    for (i in 1:length(nodeDiff)) {
      if (i == 1) {
        nat[[i]] <- list(mat[mat[, 1] == (totlen + i), 
                             2])
      }
      else {
        IN <- vector()
        P <- mat[as.numeric(mat[, 1]) <= (totlen + i), 
                 c(2, 3)]
        IN <- c(IN, P[P[, 2] == "0", 1], P[as.numeric(P[, 
                                                        2]) > (totlen + i), 1])
        nat[[i]] <- list(IN)
      }
    }
    
    # I think this is the gap to fit my code, it would make a matrix for each list in the geo_object list
    
   # for (i in 1:length(mat[, 1])) {
   #   if (mat[i, 3] == 0) {
   #     mat[i, 3] <- as.character(match(mat[i, 2], phylo$tip.label))
   #   }
   # }
   # for (i in 1:length(mat[, 1])) {
   #   if (mat[i, 3] == 0) {
   #     mat[i, 3] <- as.character(match(mat[i, 2], phylo$tip.label))
   #   }
   # }
   # if (flag == 0) {
   #   map.v <- map[match(mat[, 3], map[, 2]), 3]
   # }
   # ## I think I can insert the overlaps here, or try to call them in the third column
   # if (flag == 1) {
   #   map2 <- cbind(phylo$edge, map[, 3])
   #   map.v <- map2[match(mat[, 3], map2[, 2]), 3]
   # }
   # mat <- data.frame(mat, map.v)
   # geography.matrix <- list()
   # for (i in 1:phylo$Nnode) {
   #   var.list <- unlist(nat[[i]])
   #   len = length(var.list)
   #   int.mat <- matrix(nrow = len, ncol = len)
   #   rownames(int.mat) <- var.list
   #   colnames(int.mat) <- var.list
   #   diag(int.mat) <- 1
   #   for (j in 1:length(var.list)) {
   #     for (k in 1:length(var.list)) {
   #       if ((lower.tri(int.mat, diag = TRUE)[j, k] == 
   #            TRUE)) {
   #         int.mat[j, k] <- ifelse(mat[match(var.list[j], 
   #                                           mat[, 2]), 4] == mat[match(var.list[k], 
   #                                                                      mat[, 2]), 4], 1, 0)
   #       }
   #     }
   #   }
   #   int.mat[upper.tri(int.mat) == TRUE] <- t(int.mat)[upper.tri(t(int.mat)) == 
   #                                                       TRUE]
   #   geography.matrix[[1]] <- int.mat
   # }
    
    int.names <- rownames(vcvPhylo(phylo)) # this does not include a crown taxon representing all tips
    counter_three_letters <- 1
    for (p in (Ntip(phylo)+1):length(int.names)) {
      int.names[p] <- paste(".", THREELETTERS[counter_three_letters], 
                            sep = "")
      counter_three_letters <- counter_three_letters+1
    } # this does not g
    
    all.taxa <- unique(map[,1])
    all.sp <- list()
    all.hull <- list()

    for (p in 1:length(unique(map[,1]))) {
      #cat("iteration", p, "of", length(unique(map$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
      current.taxon <- all.taxa[[p]]
      current.data <- subset(map, map[,1] == current.taxon)
      if (nrow(current.data) > 1000) {
        current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
      } 
      points <- current.data[,c(2,3)]
      all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
      all.hull[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=1)
    }
    names(all.sp) <- unique(map[,1])
    names(all.hull) <- unique(map[,1])
    
    # make a loop that creates distributions at nodes by combining daughter distributions
    for (kk in (length(int.names)+1):(length(phylo$tip.label)+2)){
      descendantz <- Descendants(phylo, kk, type="children")
      combined.range <- gUnion(all.hull[[descendantz[1]]], all.hull[[descendantz[2]]])
      all.hull[[kk]] <- combined.range
    }
    all.hull[(Ntip(phylo)+1):(length(int.names))] <- all.hull[(Ntip(phylo)+2):(length(int.names)+1)]
    all.hull <- all.hull[1:length(int.names)]
    names(all.hull)[(Ntip(phylo)+1) : (length(all.hull))] <- int.names[(Ntip(phylo)+1):(length(int.names))]
    
    unique.combo <-  combn(int.names, m = 2)
    range_overlap <- rep(NA, ncol(unique.combo)) # make a ma
    arrange <- data.frame(species1 = unique.combo[1,], 
                          species2 = unique.combo[2,], 
                          range_overlap, stringsAsFactors=F)
    for (t in 1:length(arrange[,1])) {
      overlap <- gOverlaps(all.hull[[arrange[t,1]]], all.hull[[arrange[t,2]]])
      arrange[t,3] <- overlap
    }
    arrange[,3][arrange[,3] == "TRUE"] <- 1
    
    master.matrix <- matrix(nrow=length(int.names), ncol=length(int.names))
    master.matrix <- as.data.frame(master.matrix)
    rownames(master.matrix) <- int.names
    colnames(master.matrix) <- int.names
    #diag(master.matrix) <- 999
    for (j in 1:length(int.names)) {
      for (k in 1:length(int.names)) {
        if (j == k) {
          master.matrix[j,k] <- 1
        } else {
          master.matrix[j,k] <- filter(arrange, (species1==int.names[j] & species2==int.names[k]) |
                                   (species2==int.names[j]) & species1==int.names[k])$range_overlap
        }
      }
    }
    
    geography.matrix <- list()
    for (i in 1:phylo$Nnode) {
      var.list <- unlist(nat[[i]])
      int.mat <- master.matrix[which(rownames(master.matrix)%in%var.list),]
      int.mat <- int.mat[,which(colnames(int.mat)%in%var.list)]
      int.mat <- as.matrix(int.mat) 
      geography.matrix[[i]] <- int.mat
    }
    
   # geography.matrix <- list()
   # for (i in 1:phylo$Nnode) {
   #   var.list <- unlist(nat[[i]])
   #   len = length(var.list)
   #   int.mat <- matrix(nrow=len, ncol=len)
   #   rownames(int.mat) <- var.list
   #   colnames(int.mat) <- var.list
   #   diag(int.mat) <- 1
   #   for (j in 1:length(var.list)) {
   #     for (k in 1:length(var.list)) {
   #       if (j == k) {
   #         int.mat[j,k] <- 1
   #       } else {
   #         int.mat[j,k] <- filter(arrange, (species1==var.list[j] & species2==var.list[k]) |
   #                          (species2==var.list[j]) & species1==var.list[k])$range_overlap
   #         #int.mat[j,k] <- filter(all.matches, var.list[1]==Species1)$overlap
   #       }
   #     }
   #   }
   #   geography.matrix[[i]] <- int.mat
   # }
  }
  return(list(geography.object = geography.matrix, times = nodeDist, 
              spans = nodeDiff))
}
#test <- CreateGeoObject_SP(Anolis.data$phylo, distribution)



        
#   ## Working out the distributional overlap stuff from SD
#   int.names <- rownames(vcvPhylo(phylo)) # this does not include a crown taxon representing all tips
#   counter_three_letters <- 1
#   for (p in (Ntip(phylo)+1):length(int.names)) {
#     int.names[p] <- paste(".", THREELETTERS[counter_three_letters], 
#           sep = "")
#     counter_three_letters <- counter_three_letters+1
#   } # this does not g
#   
#   all.taxa <- unique(distribution$Name_in_Tree)
#   all.sp <- list()
#   all.hull <- list()
#   #all.poly <- list()
#   #sbbox <- make_bbox(lon = distribution$Longitude, lat = distribution$Latitude, f = .1)
#   #sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")
#   
#   for (p in 1:length(unique(distribution$Name_in_Tree))) {
#     cat("iteration", p, "of", length(unique(distribution$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
#     current.taxon <- all.taxa[[p]]
#     current.data <- subset(distribution, distribution$Name_in_Tree == current.taxon)
#     if (nrow(current.data) > 1000) {
#       current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
#     } 
#     points <- current.data[,c(2,3)]
#     all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
#     all.hull[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=1)
#     
#     # only plot the maps below if you need to check/clean the data
#     #fortified.data <- fortify(distribution.hull)
#     #pdf(paste("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Skinks_Maps/", current.taxon, ".pdf", sep=""))
#     #print(ggmap(sq_map) 
#     #  + geom_polygon(data = fortified.data, aes(x = lat, y = long, group=group, fill="tomato"), alpha=0.3)
#     #  + geom_point(data = current.data, mapping = aes(x=Longitude, y=Latitude), color="tomato3"))
#     #dev.off()
#   }
#   names(all.sp) <- unique(distribution$Name_in_Tree)
#   names(all.hull) <- unique(distribution$Name_in_Tree)
#   
#   # make a loop that creates distributions at nodes by combining daughter distributions
#   for (kk in (length(int.names)+1):(length(phylo$tip.label)+2)){
#     descendantz <- Descendants(phylo, kk, type="children")
#     combined.range <- gUnion(all.hull[[descendantz[1]]], all.hull[[descendantz[2]]])
#     all.hull[[kk]] <- combined.range
#   }
#   all.hull[(Ntip(phylo)+1):(length(int.names))] <- all.hull[(Ntip(phylo)+2):(length(int.names)+1)]
#   all.hull <- all.hull[1:length(int.names)]
#   names(all.hull)[(Ntip(phylo)+1) : (length(all.hull))] <- int.names[(Ntip(phylo)+1):(length(int.names))]
#   
#   #saveRDS(all.sp,   file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.SpatialPoints.RDS")
#   #saveRDS(all.hull, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Pygopodoids.Polygons.RDS")
#   saveRDS(polygon.list, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Pygopodoids.AllPolygons.RDS")
#   
#   
#   #range_overlap <- rep(NA, ncol(unique.combo)) # make a ma
#   #arrange <- data.frame(species1 = unique.combo[1,], 
#   #                      species2 = unique.combo[2,], 
#   #                      range_overlap, stringsAsFactors=F)
#   #for (t in 1:length(arrange[,1])) {
#   #  overlap <- gOverlaps(all.hull[[arrange[t,1]]], all.hull[[arrange[t,2]]])
#   #  arrange[t,3] <- overlap
#   #}
#   # Get unique combinations of this set:
#   inboth <- intersect(rownames(dist.mat), rownames(diff.mat)) 
#   unique.combo <-  combn(int.names, m = 2) # check ncol(unique.combo) matches length(all.TTR[[1]][,1])
#   
#   node_tip_names <- insert(int.names, (Ntip(phylo)+1), "root")
#   
#   all.matches <- NULL
#   for(p in 1:length(int.names)) {
#     phy <- phylo
#     taxon <- int.names[p] # choose a current taxon
#     
#     if (p>(Ntip(phy))) {
#       sisters <- getSisters(phy, node=(p+1), mode="number") # get the sister node/tips to this taxon
#       pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
#       pair.matrix[,1] <- taxon # who they're being compared against (taxon)
#       pair.matrix[,2] <- node_tip_names[sisters]
#       #if (sisters<=Ntip(phy)) {
#       #  pair.matrix[,2] <- rownames(dist.mat)[sisters]
#       #} else {
#       #  pair.matrix[,2] <- rownames(dist.mat)[sisters-1] # add all the comparisons
#       #}
#     } else {
#       sisters <- getSisters(phy, node=p, mode="number") # get the sister node/tips to this taxon
#       pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
#       pair.matrix[,1] <- taxon # who they're being compared against (taxon)
#       pair.matrix[,2] <- node_tip_names[sisters]
#       #if (sisters<Ntip(phy)) {
#       #  pair.matrix[,2] <- int.names[sisters]
#       #} else {
#       #  pair.matrix[,2] <- int.names[sisters] # add all the comparisons
#       #}
#     }
#     all.matches <- rbind(all.matches, pair.matrix) # keep all the pairwise matches
#   }
#   all.matches <- as.data.frame(all.matches)
#   overlap.col <- matrix(nrow=length(all.matches[,1]), ncol=1)
#   for (t in 1:length(all.matches[,1])) {
#     overlap <- gOverlaps(all.hull[[all.matches[t,1]]], all.hull[[all.matches[t,2]]])
#     overlap.col[t,1] <- overlap
#   }
#   all.matches <- cbind(all.matches, overlap.col)
#   all.matches[,3][all.matches[,3] == "TRUE"] <- 1
#   colnames(all.matches) <- c("Species1", "Species2", "overlap")
#   
#   
#   
#   
#   ### This will determine pairwise overlap for ALL tips and nodes
#   ### because of this, it takes FOREVER, so instead, I've implemented this step later
#   unique.combo <-  combn(int.names, m = 2)
#   range_overlap <- rep(NA, ncol(unique.combo)) # make a ma
#   arrange <- data.frame(species1 = unique.combo[1,], 
#                         species2 = unique.combo[2,], 
#                         range_overlap, stringsAsFactors=F)
#   for (t in 1:length(arrange[,1])) {
#     overlap <- gOverlaps(all.hull[[arrange[t,1]]], all.hull[[arrange[t,2]]])
#     arrange[t,3] <- overlap
#   }
#   arrange[,3][arrange[,3] == "TRUE"] <- 1
#   
#   