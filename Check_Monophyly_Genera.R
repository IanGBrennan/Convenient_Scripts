require(phangorn);

# check the monophyletic status of every genus in the tree `phy`
check_genus_monophyly <- function (phy) {
  # get genera. assumes form: Genus_species_whatever
  gens <- sort(unique(sub("_.*", "", phy$tip.label)));
  n <- length(gens);
  cat("Checking ", n, " genera for monophyly status.\n", sep="");
  # the number of descendant tips of the MRCA. preallocate for monotypic genera
  ndec <- rep(1, n);
  # get number of matches per genus. don't match shorter names to longer ones
  yy <- match(gsub("_.*", "_", phy$tip.label), paste0(gens, "_"));
  # the expected number of descendants == tips initially matched
  nexp <- as.numeric(table(yy));
  # get descendants for all nodes (phangorn)
  decs <- phangorn:::bip(phy); # almost instantaneous
  # skip monotypic genera
  idx <- which(nexp != 1);
  for (i in 1:length(idx)) {
  #for (i in 1:2) {
    mrca <- getMRCA(phy, which(yy == idx[i]));
    ndec[idx[i]] <- length(decs[[mrca]]);
  }
  cat("Found ", sum(ndec != nexp), " non-monophyletic genera.\n", sep="");
  return(data.frame(Genus=gens, Mono=(ndec == nexp), N.obs=ndec, N.exp=nexp, stringsAsFactors=FALSE));
}

bip <- function(x) {
  x <- reorder(x, "postorder")
  nTips <- as.integer(length(x$tip.label))
  .Call("_phangorn_bipCPP", PACKAGE = "phangorn", x$edge, nTips)
}


