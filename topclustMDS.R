inputs <- read.tree("/Users/Ian/Desktop/Species_Tree/Filtering_Trees/Full_Trees/GeneTrees>130Taxa/AllTrees.trees")
    inputs <- unroot(inputs)

    
# This function selects the best number of clusters according to the following criterion: 
# If there is a significant decline in gap from one cluster to two, then there is a single cluster; 
# otherwise, the largest increase in gap leads to the optimum number of clusters.

gapCR <- function(gaps, SEs){

      if(gaps[1] > (gaps[2] + SEs[2])){
        optimgap <- 1
      } else {
        if(any(gaps == Inf)){
                    gaps <- gaps[1:(which(gaps == Inf)[1] - 1)]
                    SEs <- SEs[1:length(gaps)]
        }
        gapdiffs <- sapply(1:(length(gaps)-1), function(i) gaps[i+1] - gaps[i])
        optimgap <- which(gapdiffs == max(gapdiffs)) + 1
        if((gaps[(optimgap - 1)] + SEs[(optimgap - 1)]) > gaps[optimgap]){
                print("The biggest difference in gaps is not significant! Returning NA.")
                return(NA)
        }
      }
      return(optimgap)
}

# This function receives a list of gene trees. It estimates the topological distance matrix. 
# MDS is then used for the given number of dimesion values and the given number of ks for each. 
# The output is a data matrix and a plot of how many dimensions were supported by each of the number 
# of dimensions.

require(cluster)
require(ape)

topclustMDS <- function(trs, mdsdim = 2, max.k, makeplot = T, criterion = "cr2", trdist = "topo"){
	
	# Create matrix of pairwise topological distances among trees
	
	topdistmat <- matrix(NA, ncol = length(trs), nrow = length(trs))
	for(i in 1:length(trs)){
	      for(j in i:length(trs)){
	      	    topdistmat[j, i] <- dist.topo(trs[[i]], trs[[j]], method = if(trdist == "topo") "PH85" else "score")
	      }
	      print(paste("finished column", i))
	}
	
	# Perform MDS representation of tree space
	
	mdsres <- cmdscale(as.dist(topdistmat), eig = T, k = mdsdim)
	if(makeplot) plot(mdsres$points)
	
	# Calculate gap statistcs for each number of clusters
	
	gapstats <- clusGap(mdsres$points, pam, max.k)
	gapstable <- gapstats$Tab
	gapstable[is.infinite(gapstable)] <- NA
        
	# The following section identifies the best-fitting number of clusters and the performs clustering usin the PAM function
	
	if(criterion == "max"){
		bestNclust <- which(gapstable[, "gap"] == max(gapstable[, "gap"], na.rm = T))
		if(length(bestNclust) > 1) bestNclust <- bestNclust[1]
	} else if(criterion == "cr2"){
	        bestNclust <- gapCR(gapstable[,3], gapstable[,4])
	}
	
	clusterdata <- pam(mdsres$points, k = bestNclust)
	
	res <- list(topodists = topdistmat, mds = mdsres, gapstats = gapstats, loci = names(data), clustering.data = clusterdata, k = bestNclust)
	
	return(res)
	
}

initial <- topclustMDS(inputs, mdsdim=2, max.k=10, makeplot = T, criterion = "cr2", trdist = "topo")

# Extract the clusters if they're of interest
cluster1 <-inputs[which(initial$clustering.data$clustering==1)]
cluster2 <-inputs[which(initial$clustering.data$clustering==2)]

# Write the trees in clusters to a separate file:
write.tree(cluster1, "/Users/Ian/Desktop/Species_Tree/Subgenera_Analysis/Varanus_genetrees_3D_CLUSTER1.trees")
write.tree(cluster2, "/Users/Ian/Desktop/Species_Tree/Subgenera_Analysis/Varanus_genetrees_3D_CLUSTER2.trees")

# Put the clusters into an object we can plot
mdsDF <- as.data.frame(initial$mds$points)
mdsDF[,3] <- as.data.frame(initial$clustering.data$clustering)
colnames(mdsDF) <- c("dim1", "dim2","cluster")

mdsDF[which(mdsDF$cluster==1),3] <- mpalette[3]
mdsDF[which(mdsDF$cluster==2),3] <- mpalette[4]



dims2 <- (ggplot(mdsDF)
          + geom_point(aes(x=dim1, y=dim2), color=mdsDF$cluster, size=3)
          + theme_classic())




all.tree <- read.tree("/Users/Ian/Desktop/Species_Tree/Filtering_Trees/Full_Trees/T474_L349.tre")
sample.list <- all.tree$tip.label

# determine how many taxa are found in all the trees provided
for (j in 1:length(inputs)){
  sample.list <- intersect(sample.list, inputs[[j]]$tip.label)
  print(paste("tree",j,":",length(sample.list), "taxa"))
}
# trim the input trees down to just the taxa found in all
new.inputs <- NULL
for (j in 1:length(inputs)){
  new.inputs[[j]] <- drop.tip(inputs[[j]], tip=setdiff(inputs[[j]]$tip.label, sample.list))
} 
class(new.inputs) <- "multiPhylo"; new.inputs <- unroot(new.inputs)

initial <- topclustMDS(new.inputs, mdsdim=2, max.k=10, makeplot = T, criterion = "cr2", trdist = "topo")


