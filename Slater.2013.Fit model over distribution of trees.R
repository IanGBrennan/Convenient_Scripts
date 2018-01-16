##### fit models over posterior distribution of trees
require(geiger);
require(fitContinuousMCMC);
source('fitContinuous_paleoModels.R');


pf <- read.table("pfile.txt", header = T); # p.file contains branch length scalars. be sure to remove first row to be able to read in file

tf <- read.nexus("tfile.nex") # tree.file

###############################################################

all.data <- read.csv("mammal_mass_posterior.csv", row.names=1); # renamed mass data to account for unedited names (i.e. Meredith et al names) in posterior distribution

all.data <- treedata(tf[[1]], all.data)$data # remove non-overlapping data

jt <- treedata(tf[[1]], all.data)$phy

all.data <- all.data[match(jt$tip.label, rownames(all.data)),] ## reorder data to be in same order as tip labels in tree 

rm(jt);

n.trees <- 100; ## fit to 100 trees

## set up matrices to hold akaike weights for all models to all datasets 

all.res <- matrix(data = NA, nrow = n.trees, ncol = 8);
extant_foss.res <- matrix(data = NA, nrow = n.trees, ncol = 8);
extant.res <- matrix(data = NA, nrow = n.trees, ncol = 8);
colnames(all.res) <- colnames(extant_foss.res) <- colnames(extant.res) <-c("bm", "trend"," acdc", "ou", "wn", "shift", "release", "randr")


trees.2.use <- sample(seq(2500, 10001), n.trees, replace = F); ## use trees after 25% burnin

rownames(all.res) <- rownames(extant_foss.res) <- rownames(extant.res) <- trees.2.use;

for(i in 1:length(trees.2.use)) {
	
	# all taxa #
	x <- trees.2.use[i]
	phy <- tf[[x]];
	phy <- treedata(phy, all.data)$phy;
	phy$edge.length <- phy$edge.length / pf$Clockrate.all.[[x]]
	weights.all <- do.model.fit(phy, all.data);
	all.res[i, ] <- weights.all;
	
	# extant and their fossil relatives #
	phy2 <- drop.tip(phy, setdiff(phy$tip.label, node.leaves(phy, mrca.of.pair(phy, "Canidae", "Tachyglossus_aculeatus"))))
	all.data2 <- treedata(phy2, all.data)$data;
	weights.2 <- do.model.fit(phy2, all.data2);
	extant_foss.res[i, ] <- weights.2;	
	
	## extant only #
	
	phy3 <- drop.tip(phy2, setdiff(phy2$tip.label, names(which(BranchingTimesFossil(phy2)==0))))
	all.data3 <- treedata(phy3, all.data2)$data;
	weights.3 <- do.model.fit(phy3, all.data3);
	extant.res[i, ] <- weights.3;
	
}

write.csv(all.res, "all_res_posterior.csv");
write.csv(extant_foss.res, "extant_foss_posterior.csv");
write.csv(extant.res, "extant.res_posterior.csv");
