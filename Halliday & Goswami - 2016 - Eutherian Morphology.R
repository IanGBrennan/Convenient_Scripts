library(ape)
library(geiger)
library(Claddis)

###############################################################
#       ****  *  *  *  *   *** ***** *  ***   *  *   ***      #
#       *     *  *  ** *  *      *   * *   *  ** *  *         #
#       ****  *  *  * **  *      *   * *   *  * **   **       #
#       *     *  *  *  *  *      *   * *   *  *  *     *      #
#       *     ****  *  *   ***   *   *  ***   *  *  ***       #
###############################################################

###  BINNING DISPARITY MATRICES ###


# For all of the bins we are interested in

BinnedDispMatrix <- function(disparitymatrix, occurrences, filename) {
  
  sapply(1:nrow(Ranges), function(x) {
    
    # Take the distance matrix of interest (whether PCO/Ancestral State/etc)
    
    disp.matrix<-disparitymatrix
    
    # For all of the edges we are interested in
    
    for (i in 353:1) {
      
      # If the edge is not in the time bin
      
      if(occurrences[i,x]==0) {
        
        # Drop that edge from the matrix
        
        disp.matrix<-disp.matrix[-i,]
        
        # But leave it in if it is in that time bin
        
      } else disp.matrix
      
    }
    
    # Save this as a matrix as part of a list of matrices
    
    write.table(disp.matrix, paste(filename,'_',x,'.txt',sep=""),
                quote=FALSE,sep='\t',col.names=TRUE, row.names=TRUE)
    
  }
  )
  
}



## ALTERNATIVE FOR DISTANCE MATRIX ##

# For all of the bins we are interested in

BinnedDispMatrix <- function(disparitymatrix, occurrences, filename) {
  
  sapply(1:nrow(Ranges), function(x) {
    
    # Take the distance matrix of interest (whether PCO/Ancestral State/etc)
    
    disp.matrix<-disparitymatrix
    
    # For all of the edges we are interested in
    
    for (i in 353:1) {
      
      # If the edge is not in the time bin
      
      if(occurrences[i,x]==0) {
        
        # Drop that edge from the matrix
        
        disp.matrix<-disp.matrix[-i,-i]
        
        # But leave it in if it is in that time bin
        
      } else disp.matrix
      
    }
    
    # Save this as a matrix as part of a list of matrices
    
    write.table(disp.matrix, paste(filename,'_',x,'.txt',sep=""),
                quote=FALSE,sep='\t',col.names=TRUE, row.names=TRUE)
    
  }
  )
  
}


### To Bootstrap the Binned Distance Matrices ###

BootstrapDistanceMatrix <- function(disparitymatrix, nrep) {
  
  # Read in the disparity matrix from the file
  
  dispmatrix<-read.table(disparitymatrix, header=T, row.names=1)
  
  # Generate a list for the matrices to go in
  
  bootlist<-list()
  
  # For as many bootstrap replicates as you want
  
  for (i in 1:nrep) {
    
    # Generate a matrix for your bootstrapped taxa
    
    bootmatrix<-matrix(nrow=nrow(dispmatrix), ncol=nrow(dispmatrix))
    
    # Randomly sample taxa (rows of the distance matrix)
    
    x<-sample(1:nrow(dispmatrix), size=nrow(dispmatrix), replace=T)
    
    # For each coordinate in the new matrix
    
    for (j in 1:nrow(dispmatrix)) {
      
      for (k in 1:nrow(dispmatrix)) {
        
        # Get the relevant distance from the original matrix
        
        bootmatrix[j,k]<-dispmatrix[x[j],x[k]]
        
      }
      
    }
    
    # Add the bootstrapped matrix to your list
    
    bootlist[[i]]<-bootmatrix
    
  }
  
  # Generate a list for the Mean Pairwise Distances
  
  MPDlist<-matrix(ncol=length(bootlist), nrow=1)
  
  # For all matrices in the bootstrapped list
  
  for (i in 1:length(bootlist)) {
    
    # Remove one half of the
    bootlist[[i]][lower.tri(bootlist[[i]], diag=T)]<-NA
    values<-na.omit(c(bootlist[[i]]))
    
    MPDlist[,i]<-mean(values)
    
  }
  
  confints<-quantile(MPDlist, c(.05, .95))
  
  output<-list(confints, MPDlist)
  
  return(output)
  
}

### To bootstrap the PCO coordinates ###

BootstrapCoordinates <- function(coordinates, nrep) {
  
  # Read in the disparity matrix from the file
  
  coords<-read.table(coordinates, header=T, row.names=1)
  
  # Generate a list for the matrices to go in
  
  bootlist<-list()
  
  
  for (i in 1:nrep) {
    
    # Generate a matrix for the coordinates to go in
    
    bootcoords<-matrix(nrow=nrow(coords), ncol=ncol(coords))
    
    
    # Randomly sample taxa (rows of the matrix)
    
    x<-sample(1:nrow(coords), size=nrow(coords), replace=T)
    
    # For each coordinate in the new coordinate
    
    for (j in 1:nrow(coords)) {
      
      for (k in 1:ncol(coords)) {
        
        # Get the relevant distance from the original matrix
        
        bootcoords[j,k]<-coords[x[j],k]
        
        
      }
      
      
    }
    
    # Add the bootstrapped matrix to your list
    
    bootlist[[i]]<-bootcoords  
  }
  
  
  
  # Generate a list for the Sums of Variances
  
  SOVlist<-matrix(ncol=length(bootlist), nrow=1)
  
  # For all matrices in the bootstrapped list
  
  for (i in 1:length(bootlist)) {
    
    # Calculate the variances of each column
    
    vars<-apply(bootlist[[i]], 2, var)
    
    # Sum those variances
    
    SOVlist[,i]<-sum(vars)
    
  }
  
  # Generate a list for the Sums of Ranges
  
  SORlist<-matrix(ncol=length(bootlist), nrow=1)
  
  # For all matrices in the bootstrapped list
  
  for (i in 1:length(bootlist)) {
    
    # Calculate the maximum of each column
    
    maximum<-apply(bootlist[[i]], 2, max)
    minimum<-apply(bootlist[[i]], 2, min)
    colranges<-maximum-minimum
    
    # Sum those ranges
    
    SORlist[,i]<-sum(colranges)
    
  }
  
  varconfints<-quantile(SOVlist, c(.05, .95))
  ranconfints<-quantile(SORlist, c(.05, .95))
  
  
  output<-list(varconfints, ranconfints)
  
  return(output)
  
}



############   DISPARITY ANALYSIS   ##############



# Ancestral State Reconstruction for trees to get node values

# A worked example, assuming that you upload a data matrix called 'matrix' and a set of topologies called 'trees'

# Create a dummy matrix to put the ancestral state matrix into (say, for DFO topology)

DFOancmatrix<-matrix

# Reconstruct Ancestral states

for (i in 1:ntrees(trees)) {
  
  DFOancmat<-AncStateEstMatrix(matrix, trees[[i]], estimate.allchars=TRUE, estimate.tips=TRUE)
  
  DFOancmatrix$matrix<-DFOancmat
  
  DFO.anc.dist.matrix<-MorphDistMatrix(DFOancmatrix)
  
  DFO.pco.data<-cmdscale(DFO.anc.dist.matrix$max.dist.matrix, k=nrow(DFO.anc.dist.matrix$max.dist.matrix) -1, add=T)$points
  
  DFOancstates<-list(DFOancstates, DFOancmat)
  
}


########### BINNING BRANCHES TO MEASURE DISPARITY  ##############


### STEP ONE - GET A MATRIX WHICH REPRESENTS THE BINS IN WHICH EACH NODE SHOULD GO ###

tree<-trees[[1]]
taxonnames<-DFOancmat

# Get the node ages

tree.node.ages<-GetNodeAges(tree)

# Create vectors to store information

branch.start<- rep(0, length(tree$edge[,1]))
branch.end<- rep(0, length(tree$edge[,1]))

# Make a matrix to put node occurrences in

nodes.in.bins<-matrix(nrow=nrow(taxonnames), ncol=nrow(Ranges))
colnames(nodes.in.bins)<-rownames(Ranges)
rownames(nodes.in.bins)<-rownames(taxonnames)
rownames(nodes.in.bins)[1:length(tree$tip.label)]<-tree$tip.label

# For each edge

for (i in 1:length(tree$edge)) {
  
  # Find the dates of either end of the edge
  
  branch.start[i]<-tree.node.ages[tree$edge[i,1]]
  branch.end[i]<-tree.node.ages[tree$edge[i,2]]
  
  
  # For each time bin
  
  for (j in 1:length(Ranges$start_time)) {
    
    # If beginning of the edge is before the end of the time bin
    
    if(branch.start[i]<Ranges$end_time[j]) {
      
      # Add a zero to the list of nodes for that bin
      nodes.in.bins[tree$edge[i,2],j]<-0 
      
      
    } else { 
      
      # If the end of the edge is after the time bin
      
      if(branch.end[i]>Ranges$start_time[j]) {
        
        # Add a zero to the list of nodes for that bin
        nodes.in.bins[tree$edge[i,2],j]<-0 
        
        
        # If neither of the above is true, the edge overlaps with the time bin			
        
      }  else {
        
        # Add a one to the list of nodes for that bin
        
        nodes.in.bins[tree$edge[i,2],j]<- 1
        
      }
    }
    
  }
  
}

# Add in the root to the oldest time bin only

nodes.in.bins[length(tree$tip.label)+1,1]<-1

for (i in 2:nrow(Ranges)) {
  nodes.in.bins[length(tree$tip.label)+1,i]<-0
}

# Order the rows of nodes.in.bins to match the matrix of interest

nodes.in.bins<-nodes.in.bins[rownames(taxonnames),,drop=FALSE]

##### USE THIS FOR BINNING A DISTANCE MATRIX #####
BinnedDistanceMatrix(taxonnames, nodes.in.bins, "Filename")

##### USE THIS FOR BINNING A PCO MATRIX #####
BinnedDispMatrix(taxonnames, nodes.in.bins, "Filename")

##### USE THIS FOR BOOTSTRAPPING A DISTANCE MATRIX #####
BootstrapDistanceMatrix(taxonnames, 1000)

##### USE THIS FOR BOOTSTRAPPING A PCO MATRIX #####
BootstrapCoordinates(taxonnames, 1000)