library(geiger)

#######################################
#### Rescale a Tree by any Variable
#######################################

tree<-read.tree("Oz.Crayfish.tre")
times<-branching.times(tree) #have a look at the branching times initially, just to check the ages
short<-rescale(tree, "depth", 49.3) #technically, this should work, though it doesn't appear the case, the function is getting masked somewhere
short<-geiger::rescale(tree, "depth", 105) #so, let's force geiger to carry out the function instead
#we used the method "depth" to rescale the deepest node to a specific value (105), this can be done with lambda, model type, etc
write.tree(short, file="Oz.Crayfish.tre") #save the tree externally

