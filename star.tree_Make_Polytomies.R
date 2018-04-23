

star.tree <- function(phy) {
  if (class(phy)=="phylo") {
    tree.string <- NULL
    for (y in 1:Ntip(phy)) {
      current.tip <- phy$tip.label[[y]]
      if (y == Ntip(phy)){
        tip.text <- paste0(current.tip, ":", max(nodeHeights(phy)))
        tree.string <- paste0(tree.string, tip.text)
        total.tree <- paste0("(", tree.string, ");")
      } else {
        tip.text <- paste0(current.tip, ":", max(nodeHeights(phy)), ",")
        tree.string <- paste0(tree.string, tip.text)
      }
    }
    phy.star <- read.tree(text=total.tree)
  } else if (class(phy)=="multiPhylo") {
    phy.star <- list()
    for (k in 1:length(phy)) {
      tree.string <- NULL
      for (y in 1:Ntip(phy[[k]])) {
        current.tip <- phy[[k]]$tip.label[[y]]
        if (y == Ntip(phy[[k]])){
          tip.text <- paste0(current.tip, ":", max(nodeHeights(phy[[k]])))
          tree.string <- paste0(tree.string, tip.text)
          total.tree <- paste0("(", tree.string, ");")
        } else {
          tip.text <- paste0(current.tip, ":", max(nodeHeights(phy[[k]])), ",")
          tree.string <- paste0(tree.string, tip.text)
        }
      }
      int.phy <- read.tree(text=total.tree)
      phy.star[[k]] <- int.phy
    }
  }

  return(phy.star)
}

# EXAMPLE:
#   normal.tree <- rtree(n=26)
#   normal.tree$tip.label <- LETTERS 
#   star.out <- star.tree(normal.tree)
#   plot(star.out)
