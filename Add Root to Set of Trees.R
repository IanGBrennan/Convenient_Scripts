## If you need to add a root for the trees, you can use this loop
## but if you write the trees, you can't read em back in
tree.list <- list()
for (i in 1:length(phymaturus.root)) {
  max.phym <- max(nodeHeights(phymaturus[[i]]))
  max.root <- max(nodeHeights(phymaturus.root[[i]]))
  diff.age <- max.root - max.phym
  tree.list[[i]] <- addroot(phymaturus[[i]], diff.age)
}
class(tree.list) <- "multiPhylo" # this is important
write.tree(tree.list, file="Phymaturus.plus.root.100.trees") # this bit, not so much, as they're useless