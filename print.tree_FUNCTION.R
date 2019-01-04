print.tree <- function (t, printflag = TRUE, relabel_nodes = FALSE, time_bp_digits = 7, 
          add_root_edge = TRUE, get_tipnames = FALSE, fossils_older_than = 0.6) 
{
  if ("node.label" %in% attributes(t)$names == FALSE) {
    rootnum = get_nodenum_structural_root(t)
    new_node_labels = paste("inNode", rootnum:(rootnum + 
                                                 t$Nnode - 1), sep = "")
    t$node.label = new_node_labels
  }
  if (relabel_nodes == TRUE) {
    rootnum = get_nodenum_structural_root(t)
    new_node_labels = paste("inNode", rootnum:(rootnum + 
                                                 t$Nnode - 1), sep = "")
    t$node.label = new_node_labels
  }
  labels = c(t$tip.label, t$node.label)
  ordered_nodenames = get_nodenums(t)
  node.types1 = rep("tip", length(t$tip.label))
  node.types2 = rep("internal", length(t$node.label))
  node.types2[1] = "root"
  node.types = c(node.types1, node.types2)
  parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, 
                                                            t$edge[, 2])
  brlen_to_parent = t$edge.length[parent_branches]
  parent_nodes = t$edge[, 1][parent_branches]
  daughter_nodes = lapply(ordered_nodenames, get_daughters, 
                          t)
  root_nodenum = get_nodenum_structural_root(t)
  tmpstr = paste("prt(t): root=", root_nodenum, "\n", sep = "")
  #prflag(tmpstr, printflag = printflag)
  levels_for_nodes = unlist(lapply(ordered_nodenames, get_level, 
                                   t))
  hts_at_end_of_branches_aka_at_nodes = t$edge.length
  hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
  h = hts_at_end_of_branches_aka_at_nodes
  times_before_present = get_max_height_tree(t) - h
  edge_ages = t$edge
  edge_ages[, 1] = h[t$edge[, 1]]
  edge_ages[, 2] = h[t$edge[, 2]]
  edge_times_bp = t$edge
  edge_times_bp[, 1] = times_before_present[t$edge[, 1]]
  edge_times_bp[, 2] = times_before_present[t$edge[, 2]]
  if (get_tipnames == TRUE) {
    list_of_clade_members_lists = rep(list(NA), length(ordered_nodenames))
    list_of_clade_members_lists[1:length(t$tip.label)] = t$tip.label
    list_of_clade_members_lists
    nontip_nodenums = (length(t$tip.label) + 1):length(ordered_nodenames)
    if (length(nontip_nodenums) > 1) {
      nontip_nodenames = ordered_nodenames[nontip_nodenums]
      nontip_cladelists = sapply(X = nontip_nodenames, 
                                 FUN = get_all_daughter_tips_of_a_node, t = t)
      nontip_cladelists
      nontip_cladelists_alphabetical = sapply(X = nontip_cladelists, 
                                              FUN = sort)
      nontip_cladelists_alphabetical
      nontip_cladelists_alphabetical_str = sapply(X = nontip_cladelists_alphabetical, 
                                                  FUN = paste, collapse = ",")
      nontip_cladelists_alphabetical_str
      list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
      list_of_clade_members_lists
    }
    else {
      nontip_nodenames = ordered_nodenames[nontip_nodenums]
      nontip_cladelists = sapply(X = nontip_nodenames, 
                                 FUN = get_all_daughter_tips_of_a_node, t = t)
      nontip_cladewords = unlist(sapply(X = nontip_cladelists, 
                                        FUN = strsplit, split = ","))
      nontip_cladelists_alphabetical = sort(nontip_cladewords)
      nontip_cladelists_alphabetical
      nontip_cladelists_alphabetical_str = paste(nontip_cladelists_alphabetical, 
                                                 collapse = ",", sep = "")
      nontip_cladelists_alphabetical_str
      list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
      list_of_clade_members_lists
    }
  }
  fossils = times_before_present > fossils_older_than
  tmpnodenums = (length(t$tip.label) + 1):(length(t$tip.label) + 
                                             t$Nnode)
  fossils[tmpnodenums] = NA
  if (get_tipnames == FALSE) {
    tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, 
                   levels_for_nodes, node.types, parent_branches, brlen_to_parent, 
                   parent_nodes, daughter_nodes, h, round(times_before_present, 
                                                          digits = time_bp_digits), fossils, labels)
    dtf = as.data.frame(tmpdtf, row.names = NULL)
    names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", 
                   "parent_br", "edge.length", "ancestor", "daughter_nds", 
                   "node_ht", "time_bp", "fossils", "label")
    dtf = unlist_dtf_cols(dtf, printflag = FALSE)
  }
  else {
    tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, 
                   levels_for_nodes, node.types, parent_branches, brlen_to_parent, 
                   parent_nodes, daughter_nodes, h, round(times_before_present, 
                                                          digits = time_bp_digits), fossils, labels, list_of_clade_members_lists)
    dtf = as.data.frame(tmpdtf, row.names = NULL)
    names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", 
                   "parent_br", "edge.length", "ancestor", "daughter_nds", 
                   "node_ht", "time_bp", "fossils", "label", "tipnames")
    dtf = unlist_dtf_cols(dtf, printflag = FALSE)
  }
  if ((add_root_edge == TRUE) && (!is.null(t$root.edge))) {
    root_row_TF = dtf$node.type == "root"
    root_edge_length = t$root.edge
    dtf$edge.length[root_row_TF] = root_edge_length
    dtf$node_ht = dtf$node_ht + root_edge_length
  }
  #prflag(dtf, printflag = printflag)
  return(dtf)
}

test <- print.tree(post.trees[[1]])
