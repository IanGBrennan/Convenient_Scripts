## This quick function pulls out the descendant tips from and edge number
#########################################################################
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
#########################################################################