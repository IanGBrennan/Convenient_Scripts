 

library(ape)

# Test it on Pyron and Burbrink
tree <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T203_Eulamprus/Constraints/PyronWBurbrink2013.tre")
ahe <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T203_Eulamprus/Constraints/Eulamprus.concat.AHE.new.names.tre")
nd4 <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T203_Eulamprus/Constraints/Eulamprus.concat.ND4.unconstrained.with.EL3882.tre")
# Let's make tree smaller for speed: drop to 400 random tips:

tree <- drop.tip(tree, tip = setdiff(tree$tip.label, sample(tree$tip.label, 400)))

  
vv <- vcv.phylo(tree)
pathsum <- diag(vv)

# give each value the correct name. This is important.
names(pathsum) <- colnames(vv)

# pathsum is now a vector of root-to-tip distances.
#   and you could compare to corresponding values for some other tree. 

# If X_UCE was a vector of root-to-tip for UCE data
#  and X_AHE was the same for AHE data, you just do 
# 
#   in_both <- intersect(names(X_UCE), names(X_AHE))
#   and you can extract values for an identical taxon set:
#    e.g., 
#     results <- data.frame(species = in_both, uce_dist = X_UCE[ in_both ], ahe_dist = X_AHE[ in_both ])
tip.dist <- data.frame(species = in_both, ahe_dist = )


#  Here I illustrate how to compute pairwise distances between taxa 
#    from two trees that potentially differ in which taxa they contain
#    You can just modify this to plot/compare root-to-tip distances 
#        instead of pairwise, using the diag(vcv.phylo(...)) trick above
#       

# To compute pairwise distances between taxa:

# And here is a matrix of pairwise distances
pairwise <- cophenetic.phylo(tree)

# if you wanted to convert this from a table to a 3 column dataframe:
ff <- as.data.frame(as.table(pairwise))


# Lets say you have a second phylogeny with different taxa and branch lengths
#   here we wil make a dummy copy of the original tree but make the branch 
#   lengths random noise

dummy <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T203_Eulamprus/Constraints/PyronWBurbrink2013.tre")
dummy$edge.length <- runif(length(dummy$edge.length))

# drop random tips
dummy <- drop.tip(dummy, tip = setdiff(dummy$tip.label, sample(dummy$tip.label, 400)))

# Get the pairwise distance matrix
pwdist_dum <- cophenetic.phylo(dummy)

# now: compare pairwise distances from these 2 trees
#   the complication is that they have different sets of taxa

# First, set of tips in common to both:
inboth <- intersect(dummy$tip.label, tree$tip.label)

# Get unique combinations of this set:
ucomb <-  combn(inboth, m = 2)

# make vectors to hold results
dist_dummy <- rep(NA, ncol(ucomb))
dist_true <-  rep(NA, ncol(ucomb))

dff <- data.frame(species1 = ucomb[1,], species2 = ucomb[2,] , dist_true, dist_dummy , stringsAsFactors=F)

# fill in the blanks....

for (ii in 1:nrow(dff)){
	
	dff$dist_true[ii]  <- pairwise[ dff$species1[ii], dff$species2[ii] ]
	dff$dist_dummy[ii] <- pwdist_dum[ dff$species1[ii], dff$species2[ii]  ]    
	
}


# plotting them... 
# Surprisingly, they should be correlated,  even with random data! 
#   this occurs due to the node density effect! 

plot(dff$dist_true, dff$dist_dummy)





ahe <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T203_Eulamprus/Constraints/Eulamprus.concat.AHE.new.names.tre")
nd4 <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T203_Eulamprus/Constraints/Eulamprus.concat.ND4.unconstrained.with.EL3882.tre")
nd4 <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T203_Eulamprus/Constraints/RAxML_bipartitions.FINAL.constrained.bsBEST.tre")
nd4 <- root(nd4, "E_egregius", resolve=T)

# First, set of tips in common to both:
inboth <- intersect(ahe$tip.label, nd4$tip.label)
drop.ahe <- setdiff(ahe$tip.label, inboth)
ahe <- drop.tip(ahe, tip=drop.ahe)
drop.nd4 <- setdiff(nd4$tip.label, inboth)
nd4 <- drop.tip(nd4, tip=drop.nd4)
plot(nd4)

# Get the pairwise distance matrices
pw.ahe <- cophenetic.phylo(ahe)
pw.nd4 <- cophenetic.phylo(nd4)

# now: compare pairwise distances from these 2 trees
#   the complication is that they have different sets of taxa

# Get unique combinations of this set:
ucomb <-  combn(inboth, m = 2)

# make vectors to hold results
dist_ahe <- rep(NA, ncol(ucomb))
dist_nd4 <- rep(NA, ncol(ucomb))

dff <- data.frame(species1 = ucomb[1,], species2 = ucomb[2,] , dist_ahe, dist_nd4, stringsAsFactors=F)

# fill in the blanks....

for (ii in 1:nrow(dff)){
  
  dff$dist_ahe[ii] <- pw.ahe[ dff$species1[ii], dff$species2[ii] ]
  dff$dist_nd4[ii] <- pw.nd4[ dff$species1[ii], dff$species2[ii] ]    
  
}

# plotting them... 
# Surprisingly, they should be correlated,  even with random data! 
#   this occurs due to the node density effect! 

plot(dff$dist_ahe, dff$dist_nd4)
factor <- (as.data.frame(dff$dist_nd4/dff$dist_ahe))
difference <- as.data.frame(dff$dist_nd4 - dff$dist_ahe)
dff <- cbind(dff, factor[,1])
dff <- cbind(dff, difference)
densityplot(dff[,5])
test<-densityplot(dff[,6])

dff$group <- "test"

p <- ggplot(dff, aes(group, `factor[, 1]`))
p + geom_violin() + geom_jitter(height = 0, width = 0.1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))


#######################################################################################
# Interlude: the base 'plot' and 'abline' functions are alright, but we want to 
## make it (1) prettier, and (2) include the information from our linear regression
### into the plot, so that we know what our results were. Use custom 'ggplotRegression'
### if you want to change the saturation use 'alpha'
ggplotRegression <- function (fit) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(alpha=0.25, color="red") + # change to 0.25 and "red" for time plots
    stat_smooth(method = "lm", col = "black") + # change to "black" for time plots
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
#######################################################################################
fit <- lm(dist_nd4 ~ dist_ahe, data=dff) # change this according to the parameter you simulated
plot.fit <- (ggplotRegression(fit))
plot.fit + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))










####### Root-Tip Distances

ahe_vv <- vcv.phylo(ahe)
nd4_vv <- vcv.phylo(nd4)

ahe_path <- diag(ahe_vv)
nd4_path <- diag(nd4_vv)

# give each value the correct name. This is important.
names(ahe_path) <- colnames(ahe_vv)
names(nd4_path) <- colnames(nd4_vv)

# pathsum is now a vector of root-to-tip distances.
#   and you could compare to corresponding values for some other tree. 

# If X_UCE was a vector of root-to-tip for UCE data
#  and X_AHE was the same for AHE data, you just do 
# 
#   in_both <- intersect(names(X_UCE), names(X_AHE))
#   and you can extract values for an identical taxon set:
#    e.g., 
#     results <- data.frame(species = in_both, uce_dist = X_UCE[ in_both ], ahe_dist = X_AHE[ in_both ])
tip.dist <- data.frame(species = inboth, ahe_dist = ahe_path[inboth], nd4_dist = nd4_path[inboth])
plot(tip.dist$ahe_dist, tip.dist$nd4_dist)
factor.dist <- (as.data.frame(tip.dist$nd4_dist/tip.dist$ahe_dist))
tip.dist <- cbind(tip.dist, factor.dist[,1])
densityplot(tip.dist[,4])
