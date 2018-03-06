tree<-pbtree(n=100,scale=1)
x<-fastBM(phytools:::ebTree(tree,-4),internal=T)
a<-x[1:tree$Nnode+Ntip(tree)]
x<-x[1:Ntip(tree)]
phenogram(tree,c(x,a),ftype="off",col="lightblue")

fitEB<-anc.ML(tree,x,model="EB")
fitEB$sig2

plot(tree)
tree <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Agamids.tre")
simmap.tree<-make.era.map(tree,c(0,max(nodeHeights(tree)-10)))
plotSimmap(simmap.tree)
test <- map.to.singleton(simmap.tree)
plotTree.singletons(test)


test <- chainsaw2(tree, 20, return_pieces = F)
test
plot(test)
for (i in 1:Ntip(test)) {
  if (test$tip.label[[1]])
}

chop.tree <- function (tr, timepoint = 10, return_pieces = TRUE) 
{
  tr_table = prt(tr, printflag = FALSE, get_tipnames = FALSE)
  tr_table
  TF_exists_more_recently_than_10mya = tr_table$time_bp < timepoint
  labels_for_tips_existing_more_recently_than_10mya = tr_table$label[TF_exists_more_recently_than_10mya == 
                                                                       TRUE]
  edge_times_bp = get_edge_times_before_present(tr)
  edges_start_earlier_than_10mya = edge_times_bp[, 1] > timepoint
  edges_end_later_than_10mya = edge_times_bp[, 2] <= timepoint
  edges_to_chainsaw = edges_start_earlier_than_10mya + edges_end_later_than_10mya == 
    2
  nodes_to_chainsaw = tr$edge[, 2][edges_to_chainsaw]
  numtips = length(tr$tip.label)
  tree_to_chainsaw = tr
  if (return_pieces == TRUE) {
    return_pieces_list = as.list(rep(NA, length(nodes_to_chainsaw)))
    return_pieces_basenames = as.list(rep(NA, length(nodes_to_chainsaw)))
    chopTable = NULL
  }
  chainsaw_table = NULL
  for (i in 1:length(nodes_to_chainsaw)) {
    if (nodes_to_chainsaw[i] <= numtips) {
      if (return_pieces == TRUE) {
        return_pieces_list[[i]] = timepoint
        tmp_tipname = tr$tip.label[nodes_to_chainsaw[i]]
        return_pieces_basenames[[i]] = tmp_tipname
      }
    }
    #else{
    #  nodes_to_chainsaw[i] <- paste("Node", nodes_to_chainsaw[i], sep=".")
    #}
    
    else {
      tmp_subtree = extract.clade(tr, nodes_to_chainsaw[i])
      branchlength_below_subtree_LCA_node = timepoint - 
        get_max_height_tree(tmp_subtree)
      tmp_subtree$root.edge = branchlength_below_subtree_LCA_node
      if (return_pieces == TRUE) {
        return_pieces_list[[i]] = tmp_subtree
        tmp_labels_merge = paste(tmp_subtree$tip.label, 
                                 collapse = ",", sep = "")
        tmp_labels_split = strsplit(tmp_labels_merge, 
                                    split = ",")[[1]]
        new_labels = sort(tmp_labels_split)
        basename_after_cutting = paste(new_labels, collapse = ",", 
                                       sep = "")
        return_pieces_basenames[[i]] = basename_after_cutting
      }
      tmp_number_of_tips = length(tmp_subtree$tip.label)
      numtips_to_drop = tmp_number_of_tips - 1
      tmp_labels = tmp_subtree$tip.label
      labels_to_drop = tmp_labels[1:numtips_to_drop]
      ordered_labels_to_make_into_new_name = sort(tmp_labels)
      #name_new_tip = paste(ordered_labels_to_make_into_new_name, 
      #                     collapse = ",", sep = "")
      name_new_tip = paste("Node", nodes_to_chainsaw[i], sep=".")
      label_kept_num = length(tmp_labels)
      label_kept = tmp_labels[label_kept_num]
      new_label = name_new_tip
      tree_to_chainsaw$tip.label[tree_to_chainsaw$tip.label == 
                                   label_kept] = new_label
      tree_to_chainsaw = drop.tip(tree_to_chainsaw, labels_to_drop)
    }
  }
  tree_to_chainsaw_table = prt(tree_to_chainsaw, printflag = FALSE)
  tree_to_chainsaw_table_tips_TF_time_bp_LT_10my = tree_to_chainsaw_table$time_bp < 
    timepoint
  tmp_edge_lengths = tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
  times_bp_for_edges_to_chainsaw = tree_to_chainsaw_table$time_bp[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
  adjustment = times_bp_for_edges_to_chainsaw - timepoint
  revised_tmp_edge_lengths = tmp_edge_lengths + adjustment
  tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my] = revised_tmp_edge_lengths
  ordered_nodenames = get_nodenums(tree_to_chainsaw)
  parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, 
                                                            tree_to_chainsaw$edge[, 2])
  NA_false = is.not.na(tree_to_chainsaw_table$edge.length)
  tree_to_chainsaw$edge.length[parent_branches[NA_false]] = tree_to_chainsaw_table$edge.length[NA_false]
  if (return_pieces == TRUE) {
    chainsaw_result = NULL
    chainsaw_result$tree_to_chainsaw = tree_to_chainsaw
    chainsaw_result$return_pieces_list = return_pieces_list
    chainsaw_result$return_pieces_basenames = return_pieces_basenames
    class(chainsaw_result) = "chainsaw_result"
    return(chainsaw_result)
  }
  else {
    return(tree_to_chainsaw)
  }
}

test <- chop.tree(tree, 10, return_pieces = F)
plot(test)
scaled <- rescale(test, "EB", -0.5)
plot(scaled)


test <- branching.times(tree)
densityplot(test)
mean(test)
mode(test)
