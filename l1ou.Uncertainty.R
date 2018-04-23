library(geiger)
library(phytools)
library(ggbiplot);library(ggplot2);library(ggtree); library(ggridges)
library(l1ou)
library(Rmisc)
library(wesanderson); library(ggthemes)

# load your trees
#trees = read.tree("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Final_Trees/Marsupials.ASTRAL.RAxML.trees") #pick out our set of posterior trees

# load your data, with headers for columns (if applicable), and rownames as the taxon names
#trait.data <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/Marsupials.RAW.csv", header=T, row.names=1)
  # there's an interesting thing to take into consideration here, if you're doing a phyloPCA or phylo-corrected regression
  # then your data will depend on the tree, and you'll need to make a loop to correct the data for each tree
  # if you're not doing some sort of phylo correction, you don't need to do this, just bear in mind

########################################################
# Now we can start looking at shifts in size and shape
## using 'l1ou' we'll estimate the position of shifts
## then attempt to identify instances of convergence.
# First we'll just do it with one tree to practice
########################################################

## We'll need this function to pull out the descendant tips and edge numbers
############################################################################
getDescendants.edges<-function(tree,edge,curr=NULL){
  names <- NULL
  if(is.null(curr)) curr<-vector()
  node.below <- tree$edge[edge,2]
  if(node.below <= Ntip(tree)) {
    input <- tree$tip.label[[node.below]]
    names <- append(names, input)
  }
  else {
    daughters<-tree$edge[which(tree$edge[,1]==node.below),2]
    curr<-c(curr,daughters)
    z<-which(daughters<=length(tree$tip))
    if(length(z)==2) for(i in 1:length(z)) {
      input <- tree$tip.label[[curr[[i]]]]
      names <- append(names, input)
    }
    if(length(z)==1) {
      target <- daughters[[z]]
      input <- tree$tip.label[[target]]
      names <- append(names, input)
    }
    w<-which(daughters>=length(tree$tip))
    if(length(w)>0) for(i in 1:length(w)) 
      curr<-getDescendants(tree,daughters[w[1]],curr)
    curr<-unique(curr)
    curr<-subset(curr, curr<=Ntip(tree))
    for (q in 1:length(curr)) {
      input <- tree$tip.label[[curr[[q]]]]
      names <- append(names, input)
    }
  }
  names <- unique(names)
  return(names)
}
############################################################################


## This function wraps up the l1ou analyses across multiple trees.
  ## Arguments
      # phy - a multiphylo object (nexus or newick)
      # traits - a data frame where rownames are taxon names, and columns (with colnames) are the traits
      # n.iter - the number of iterations (trees) to run the analysis over, cannot exceed the number of trees
                #provided in the phy object
      # estimate.convergence - TRUE or FALSE, should we estimate shifts, then determine convergent regimes?
  ## Value
      # number.shifts - data frame of the number of shifts per tree
      # shift.position.list - vector of all shift positions across all analyzed trees
      # shift.position.by.tree - list of vectors of shift positions (stored by tree, $shift.position.list[[1]] = results from tree[[1]])
      # l1ou.res - list of l1ou output values ($l1ou.res[[1]] = results from tree[[1]])

estimate.uncertainty <- function(phy, traits, n.iter=10, estimate.convergence=TRUE, output.dir=NULL) {
  no.shifts <- NULL
  shift.positions.by.tree <- list()
  shift.positions.list <- NULL
  l1ou.res <- NULL
  
  #if (n.iter > length(phy)) {
  #  print("sorry, number of iterations cannot exceed number of provided trees")
  #  stop()
  #} else 
  if (estimate.convergence==TRUE) {
    for (i in 1:n.iter) {
      cat("iteration", i, "of", n.iter, "\n") #keep track of what tree/loop# we're on
      
      #### Adjust the data and tree to fit (order matters in l1ou!)
      data <- adjust_data(phy[[i]], traits[[i]]) # specify the data columns of interest
      
      #### Estimate the number and position of shifts a priori 
      fit <- estimate_shift_configuration(data$tree, data$Y, nCores=8, quietly=F, criterion="pBIC") # if you want to do only a single trait, 'data$Y[,x]'
      shift.fit <- estimate_convergent_regimes(fit, nCores=8, criterion="pBIC")
      
      l1ou.res[[i]] <- shift.fit
      #plot(shift.fit, cex=0.8)
      # if you get 'figure margins' error, do: par(mar=c(1,1,1,1))
      
      #### Create a data frame to hold the tree # and the # of shifts inferred
      shifts.frame <- NULL
      shifts.frame <- as.data.frame(t(c(i, shift.fit$nShifts)))
      no.shifts <- rbind(no.shifts, shifts.frame)
      
      # have to pull the shift positions (edges)
      shift.edges <- shift.fit$shift.configuration
      # match it to tips
      all.shifted.tips <- NULL
      if (length(shift.edges) == 0) {
        all.shifted.tips <- "no shifts"
      } else for (t in 1:length(shift.edges)) {
        names <- getDescendants.edges(data$tree, shift.edges[[t]])
        all.shifted.tips <- append(all.shifted.tips, names)
      }
      shift.positions.list <- append(shift.positions.list, all.shifted.tips)
      shift.positions.by.tree[[i]] <- all.shifted.tips
      
      if (!is.null(output.dir)) {
        dir.work <- getwd()
        lab.no.shifts <- paste0(dir.work, "/", output.dir, "/num.shifts.RDS")
        lab.shift.positions.list <- paste0(dir.work, "/", output.dir, "/list.shift.positions.RDS")
        lab.shift.positions.by.tree <- paste0(dir.work, "/", output.dir, "/shift.positions.by.tree.RDS")
        lab.l1ou.res <- paste0(dir.work, "/", output.dir, "/l1ou.Results.RDS")
        
        saveRDS(no.shifts,               file = lab.no.shifts)
        saveRDS(shift.positions.list,    file = lab.shift.positions.list)
        saveRDS(shift.positions.by.tree, file = lab.shift.positions.by.tree)
        saveRDS(l1ou.res,                file = lab.l1ou.res)
      }
    }
  } else if (estimate.convergence==FALSE){
    for (i in 1:n.iter) {
      cat("iteration", i, "of", n.iter, "\n") #keep track of what tree/loop# we're on
      
      #### Adjust the data and tree to fit (order matters in l1ou!)
      data <- adjust_data(phy[[i]], traits[[i]]) # specify the data columns of interest
      
      #### Estimate the number and position of shifts a priori 
      shift.fit <- estimate_shift_configuration(data$tree, data$Y, nCores=8, quietly=F, criterion="pBIC") # if you want to do only a single trait, 'data$Y[,x]'
      l1ou.res[[i]] <- shift.fit
      #plot(shift.fit, cex=0.8)
      # if you get 'figure margins' error, do: par(mar=c(1,1,1,1))
      
      #### Create a data frame to hold the tree # and the # of shifts inferred
      shifts.frame <- NULL
      shifts.frame <- as.data.frame(t(c(i, shift.fit$nShifts)))
      no.shifts <- rbind(no.shifts, shifts.frame)
      
      # have to pull the shift positions (edges)
      shift.edges <- shift.fit$shift.configuration
      # match it to tips
      all.shifted.tips <- NULL
      if (length(shift.edges) == 0) {
        all.shifted.tips <- "no shifts"
      } else for (t in 1:length(shift.edges)) {
        names <- getDescendants.edges(data$tree, shift.edges[[t]])
        all.shifted.tips <- append(all.shifted.tips, names)
      }
      shift.positions.list <- append(shift.positions.list, all.shifted.tips)
      shift.positions.by.tree[[i]] <- all.shifted.tips
      
      if (!is.null(output.dir)) {
        dir.work <- getwd()
        lab.no.shifts <- paste0(dir.work, "/", output.dir, "/num.shifts.RDS")
        lab.shift.positions.list <- paste0(dir.work, "/", output.dir, "/list.shift.positions.RDS")
        lab.shift.positions.by.tree <- paste0(dir.work, "/", output.dir, "/shift.positions.by.tree.RDS")
        lab.l1ou.res <- paste0(dir.work, "/", output.dir, "/l1ou.Results.RDS")
        
        saveRDS(no.shifts,               file = lab.no.shifts)
        saveRDS(shift.positions.list,    file = lab.shift.positions.list)
        saveRDS(shift.positions.by.tree, file = lab.shift.positions.by.tree)
        saveRDS(l1ou.res,                file = lab.l1ou.res)
      }
    }
  } 
  colnames(no.shifts) <- c("tree.no", "n.shifts")
  result <- list()
  result$number.shifts <- no.shifts
  result$shift.positions.list <- shift.positions.list
  result$shift.positions.by.tree <- shift.positions.by.tree
  result$l1ou.res <- l1ou.res
  return(result)
}

#example:
#test <- estimate.uncertainty(trees, trait.data, n.iter=10, estimate.convergence=T)


## This next function will process and plot the output from our analysis above
  ## Arguments
      # estimate - the output from the 'estimate.uncertainty' function
      # tree.num - the number of the tree that you want to plot the final shift set on
  ## Value
      # all the same outputs from a single l1ou run (Y, tree, shift.configuration, etc.)

process.uncertainty <- function(estimate, tree.num) {
  chosen.tree <- estimate$l1ou.res[[tree.num]]$tree #designate tree to plot on
  res.counts <- table(estimate$shift.positions.list) # make a table of the shift frequencies and placements
  shifted.tips <- as.data.frame(res.counts) # turn it into a data frame
  colnames(shifted.tips) <- c("shift.positions.list", "Freq")
  all.node.numbers <- as.data.frame(chosen.tree$tip.label) # get all the nodes of the tree and the numbers, the tree must match the one you want to plot!
  all.node.numbers[,"tip.no"] <- rownames(all.node.numbers); colnames(all.node.numbers) <- c("tip.name", "tip.no") # make a column that shows the tip number
  target.numbers <-  all.node.numbers[all.node.numbers$tip.name %in% shifted.tips$shift.positions.list,] # subset all the tips, so show just the node numbers of the shifted tips
  target.numbers[,2] <- as.numeric(target.numbers[,2]) # make it numeric
  target.numbers <- target.numbers[order(match(target.numbers$tip.name, shifted.tips$shift.positions.list)),] # match the new frame to the output (shift) frame
  target.numbers[,"shift.freq"] <- cbind(shifted.tips$Freq) # add on the frequencies of shifts
  
  chosen.shift <- estimate$l1ou.res[[tree.num]] # designate which shift set you want to plot the adjustments on 
  tree <- chosen.shift$tree
  
  # transform the shifted branches using this loop
  for (i in 1:length(target.numbers[,2])){
    rownames(tree$edge) <- c(1:length(tree$edge[,1])) # give the tree edge frame rownames
    target.almost <- subset(tree$edge, tree$edge[,2]==target.numbers[,2][i]) # pull out the ancestor and descendant nodes of the target edge
    interim.target <- subset(target.numbers, target.numbers[,2]==target.numbers[,2][i]) # subset the data frame to just the current tip of interest (descendant node)
    target.edge <- as.numeric(rownames(target.almost)) # get the number of the target edge
    tree$edge.length[[target.edge]] <- tree$edge.length[[target.edge]]+(0.01*interim.target[,3]) # add the desired length to the branch, per shift (here, 0.01)
  }
  chosen.shift$tree <- tree # set the tree in the 'shift.fit' object to our rescaled tree
  plot(tree, show.tip.label=F)
  return(chosen.shift)
}

#example:
#test.out <- process.uncertainty(test, 1)
#plot(test.out)

# trying to find the consistency of convergent shifts
#ouchie <- matrix(nrow=length(test$l1ou.res[[1]]$shift.configuration), ncol=3)
#prac <- as.matrix(test$l1ou.res[[1]]$shift.configuration)
#ouchie[,1] <- as.numeric(rownames(prac)); ouchie[,2] <- prac[,1]
#ouchie <- as.data.frame(ouchie)
#colnames(ouchie) <- c("regime", "tip.no", "tip.name")
#for (k in 1:length(ouchie[,3])) {
#  current.tip <- ouchie[k,2]
#  name <- test$l1ou.res[[1]]$tree$tip.label[current.tip] # test = "estimate"
#  
#}


