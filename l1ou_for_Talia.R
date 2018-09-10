library(geiger)
library(phytools)
#library(ggbiplot);library(ggplot2);library(ggtree); library(ggridges)
library(l1ou)
library(Rmisc)
#library(wesanderson); library(ggthemes)

source("/PATH_TO/l1ou.Uncertainty.R")
source("/PATH_TO/Get.Descendant.Edges.R")

# load your trees
trees = read.tree("/PATH_TO/trees") #pick out our set of posterior trees

# load your data, with headers for columns (if applicable), and rownames as the taxon names
trait.data <- read.csv("/PATH_TO/data.csv", header=T, row.names=1)
# there's an interesting thing to take into consideration here, if you're doing a phyloPCA or phylo-corrected regression
# then your data will depend on the tree, and you'll need to make a loop to correct the data for each tree
# if you're not doing some sort of phylo correction, you don't need to do this, just bear in mind

# I assume you can read up on the l1ou method and package on your own, so just the code bits:

########################################################
# Now we can start looking at shifts in size and shape
## using 'l1ou' we'll estimate the position of shifts
## then attempt to identify instances of convergence.
# First we'll just do it with one tree to practice
########################################################
    #### Adjust the data and tree to fit (order matters in l1ou!)
    data <- adjust_data(tree, trait.data) # you can specify the data columns of interest with $, if phy object is multiphylo, designate a single tree with [[1]]
    
    #### Estimate the number and position of shifts a priori 
    shift.fit <- estimate_shift_configuration(data$tree, data$Y, nCores=8, quietly=F, criterion="pBIC") # if you want to do only a single trait, 'data$Y[,x]'
    plot(shift.fit) # if you get 'figure margins' error, do: par(mar=c(1,1,1,1))
    
    #### Determine if there's convergence between shifted regimes
    conv.fit <- estimate_convergent_regimes(shift.fit, nCores=8, criterion="pBIC")
    plot(conv.fit)



# If you want to account for uncertainty in phylo reconstruction, try iterating across multiple trees
    # this script runs a model fit for each tree, then adds a fractional unit of length to each branch
    # that exhibits a shift in the posterior. More shifts = more consistent support for a shift.
#####################################################################################################    
    #### Estimate the number and position of shifts a priori 
    output.est <- estimate.uncertainty(trees, trait.data, n.iter=2, estimate.convergence=TRUE)
    # remember 'trees' has to be multi.phylo, trait.data has to be a list of data frames
    
    #### Process the output from the uncertainty, plot results on the tree you prefer (MCC?)
    output.post <- process.uncertainty(output.est, tree.num = 1)
    # look at 'output.post' and it will give you the shift/convergence config for your chosen tree
    
    plot(output.post)
    # if you get 'figure margins' error, do: par(mar=c(1,1,1,1))


