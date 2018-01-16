#################################
# PATTERNS OF RATE HETEROGENEITY#
#################################
library(phytools)
library(geiger)
library(phangorn)
library(arbutus)
library(phylobase)

## functions to estimate prevalence and patterns of rate-heterogeneity: (i) clade events, (ii) rate shifts on internal branches, not passed to descendants, and (iii) rate-changes on isolated, terminal branches

# functions take: the original, input phylogeny (phy_original)
                 # a vector of tip-data (traits)
                 # the (median) scaled tree from variable-rates models (phy_rescaled)
                 # a threshold at which a rate change (proportion between the branch length on the scaled tree and the length of the identical branch in the input phylogeny)is counted a rate-shift

# read in your data
tree<-read.nexus('BT.Pygopodoidea.tre')
data<-read.csv("BT.Pygopodoidea.logSVL.csv", row.names=1, header=F)
rescaledtree<-read.nexus('BT.Pygopodoidea.logSVL.con.tre')

# define fxn 'gener': determines first X generations of node 'nodu' (here maxim = number of generations)
####################################################################################
gener <- function(tree, nodu, start, maxim) #for now, it will also include the node about which you are asking, but it doesn't matter
{ rownames(tree$edge)<- 1:length(tree$edge.length)
not_allowed<-getDescendants(tree,nodu)

lngth=2
if (length(not_allowed)<lngth)
  lngth=length(not_allowed)
returned_value= nodu # a list of not_allowed nodes i.e. first descendants
if (start < maxim) #get only the first maxim generations
  for(i in 1:lngth)
    if (not_allowed[i] != nodu) #so it won't put the same node many times
      returned_value = c(returned_value,gener(tree, not_allowed[i],start+1,maxim )) #get the children of the children 
return (returned_value)
}
####################################################################################


# define the function to estimate rate - increases (drop-down)
####################################################################################
ht_pattern_increase<-function(phy_original, traits, phy_rescaled, threshold)
{#get rates
  rates<- phy_rescaled$edge.length/phy_original$edge.length
  rates_tree<- phy_rescaled
  rates_tree$edge.length<-rates #compute br_len for more elegant code
  #return: n_clade, n_internal, prop_terminal i.e. number of clade shifts, number of internal shifts not passed to descendants, and proportion of terminal branches with shifts out of the total number of terminal branches
  # which_clade, which_internal, which_term - what is the parental node of the clade involved in the clade event & which are the internal and terminal branches that have shifts
  n_clade<-0
  which_clade<-numeric() # i.e. father node of clade
  n_internal<-0
  which_internal <- numeric() # i.e. which internal branches
  prop_terminal<-0
  which_term<-numeric() # i.e. which terminal branches
  
  #(1) clade event
  # for each non-terminal node -> search descendant clade
  # if found one, cannot search through its descendants i.e. don't double count clades within clades (the ifelse if part of the code)
  for(jj in (length(rates_tree$tip.label)+1): (rates_tree$Nnode + length(rates_tree$tip.label)) )
    # get node descendants, must be >= 4
    if(length(Descendants(rates_tree,jj,"all"))>=4)
      if( ifelse (length(which_clade)>0, !jj %in%unlist(Descendants(rates_tree,as.numeric(which_clade),"all")), TRUE))
        # if there are not nodes in which_clade, returns TRUE
        # if there are nodes in which_clade, return: false if node is in descendants of counted_clade & true if node is not in descendants of counted_clade
        if( sum(rates_tree$edge.length[match(Descendants(rates_tree,jj,"all"),rates_tree$edge[,2])]>threshold) / length(rates_tree$edge.length[match(Descendants(rates_tree,jj,"all"),rates_tree$edge[,2])])>0.75) # 75% of clade members >threshold
        {n_clade<-n_clade+1 # count clade event
         which_clade<-c(which_clade, jj)}
  
  # (2) rate changes at tips
  # get terminal branches
  rownames(rates_tree$edge)<- 1:length(rates_tree$edge.length)
  t_br<-as.numeric(rownames(subset(rates_tree$edge, rates_tree$edge[,2]<=length(rates_tree$tip.label))))
  
  # proportion of terminal branches > threshold, BUT condition: must not be part of the clade_event(s)
  term_clade<-numeric() # which terminal branches are part of the clade event(s)?
  if(length(which_clade)>0)
    for(jj in 1: length(which_clade)) 
      term_clade<-c(term_clade,as.numeric(names(rates_tree$edge[,2][rates_tree$edge[,2]%in%Descendants(rates_tree,which_clade[[jj]],"all")]))[as.numeric(names(rates_tree$edge[,2][rates_tree$edge[,2]%in%Descendants(rates_tree,which_clade[[jj]],"all")]))%in%t_br ])
  
  # prop_terminal =  terminal branches that are >threshold and not part of the clade_event(s) / total number of terminal branches
  prop_terminal<- sum(rates_tree$edge.length[t_br]>threshold & !(t_br%in%term_clade))  / length(t_br)
  which_term<- t_br[rates_tree$edge.length[t_br]>threshold & !(t_br%in%term_clade)]
  
  #(3) internal single-lineage bursts
  # internal branches > threshold & not part of the clade event - cannot be part of it anyway because then the shift would be passed to descendants
  # so condition is that int_branch > threshold and its first 3 generations cannot be bigger than threshold (i.e. shift not passed)
  # the branch can have only 2 descendants; if you want that descendants form a clade i.e. more than 4 - just code extra condition
  int_br<-(1:length(rates_tree$edge.length))[-t_br]
  for(jj in 1: length(int_br))
    if( rates_tree$edge.length[int_br[[jj]]]>threshold ) 
      if( all ( rates_tree$edge.length[tail(as.numeric(names( gener(rates_tree, rates_tree$edge[int_br[[jj]],2], 0, 3))),-1)]<threshold) )  #i.e. a high rate is not passed to descendants 
        # gener(rates_tree, rates_tree$edge[int_br[[jj]],2], 0, 3) includes the end-node of int_br[[jj]] -> should be excluded, so used 'tail'
      {n_internal<-n_internal+1  #count branch
       which_internal <- c(which_internal, int_br[[jj]])}
  
  ret<- list(n_internal, n_clade, prop_terminal, which_internal, which_clade, which_term)
  names(ret) <- c("n_internal", "n_clade", "prop_term", "which_int", "which_clade", "which_term")
  return(ret)
}
####################################################################################

# test out the rate increase fxn
increase<-ht_pattern_increase(tree, data, rescaledtree, threshold=2)
Descendants(tree, node=349)

# define the function to estimate rate-decreases (drop-down)
# takes same arguments, essentially same function but this time rates < threshold 
####################################################################################
ht_pattern_decrease<-function(phy_original, traits, phy_rescaled, threshold)
{#get rates
  rates<- phy_rescaled$edge.length/phy_original$edge.length
  rates_tree<- phy_rescaled
  rates_tree$edge.length<-rates #compute br_len would be more elegant
  n_clade<-0
  which_clade<-numeric() # i.e. father node of clade
  n_internal<-0
  which_internal <- numeric() # i.e. which internal branches
  prop_terminal<-0
  which_term<-numeric() # i.e. which terminal branches
  
  #(1) clade
  for(jj in (length(rates_tree$tip.label)+1): (rates_tree$Nnode + length(rates_tree$tip.label)) )
    if(length(Descendants(rates_tree,jj,"all"))>=4)
      if( ifelse (length(which_clade)>0, !jj %in%unlist(Descendants(rates_tree,as.numeric(which_clade),"all")), TRUE))
        if( sum(rates_tree$edge.length[match(Descendants(rates_tree,jj,"all"),rates_tree$edge[,2])]<threshold) / length(rates_tree$edge.length[match(Descendants(rates_tree,jj,"all"),rates_tree$edge[,2])])>0.75 ) # 75% of clade members <threshold
        {n_clade<-n_clade+1 # count clade event
         which_clade<-c(which_clade, jj)}
  
  # (2) terminal shifts
  rownames(rates_tree$edge)<- 1:length(rates_tree$edge.length)
  t_br<-as.numeric(rownames(subset(rates_tree$edge, rates_tree$edge[,2]<=length(rates_tree$tip.label))))
  
  term_clade<-numeric() # which terminal branches are part of the clade event(s)?
  if(length(which_clade)>0)
    for(jj in 1: length(which_clade)) 
      term_clade<-c(term_clade,as.numeric(names(rates_tree$edge[,2][rates_tree$edge[,2]%in%Descendants(rates_tree,which_clade[[jj]],"all")]))[as.numeric(names(rates_tree$edge[,2][rates_tree$edge[,2]%in%Descendants(rates_tree,which_clade[[jj]],"all")]))%in%t_br ])
  
  prop_terminal<- sum(rates_tree$edge.length[t_br]<threshold & !(t_br%in%term_clade))  / length(t_br)
  which_term<- t_br[rates_tree$edge.length[t_br]<threshold & !(t_br%in%term_clade)]
  
  #(3) internal shifts
  int_br<-(1:length(rates_tree$edge.length))[-t_br]
  for(jj in 1: length(int_br))
    if( rates_tree$edge.length[int_br[[jj]]]<threshold ) 
      if( all ( rates_tree$edge.length[tail(as.numeric(names( gener(rates_tree, rates_tree$edge[int_br[[jj]],2], 0, 3))),-1)]>threshold) ) 
      {n_internal<-n_internal+1
       which_internal <- c(which_internal, int_br[[jj]])}
  
  ret<- list(n_internal, n_clade, prop_terminal, which_internal, which_clade, which_term)
  names(ret) <- c("n_internal", "n_clade", "prop_term", "which_int", "which_clade", "which_term")
  return(ret) }
####################################################################################

#test out the rate decrease fxn
decrease<-ht_pattern_increase(tree, data, rescaledtree, threshold=0.5)

# plot the tree with the associated labels to determine where the shifts are
plotTree(tree); nodelabels(frame="circle", cex=0.2) 
plotTree(tree); tiplabels(frame="circle", cex=0.2)
plotTree(tree); edgelabels(frame="circle", cex=0.2)
