# The function requires a few things:
## phy: your tree file in phylo format
## cf.file: the path to your concordance factor file (e.g. "~/Desktop/cf.stat")
## terminal.color: terminal branches don't have CFs, so how do you want them colored? Default is the color set for CF=100
## col.palette: the color palette from RColorBrewer that you'd like. Default is "RdYlBu"
## legend: TRUE/FALSE, do you want a legend of the gCF colors?

color.by.CF <- function(phy, cf.file, terminal.color=NULL, col.palette="RdYlBu", legend=TRUE, value.check=FALSE){
  require(RColorBrewer); require(dplyr)
  
  # set the plot layout
  if(legend==T){layout(matrix(1:2,ncol=2), width = c(4,1),height = c(1,1))}
  
  # read in the gCF file
  cf.data = read.delim(cf.file, header=T, comment.char="#")
  
  # rename the first column from 'ID' to 'child', 
  # the naming convention for branches in IQTREE is not the same for edges in APE
  # but it is the same as the number of the child node below an edge in APE, so we can use that
  names(cf.data)[1] = "child"
  # add root node, and apply an arbitrary value, this may just be funny behavior because of 
  cf.data[1,"gCF"] <- 100; # arbitrary CF value for 2nd node
  cf.data <- rbind(cf.data, c(min(cf.data$child)-1, 100, NA,NA,NA,NA,NA)) # arbitrary CF value for root node
  
  # create a data frame of all the edges in the tree, listed by (parent, child) nodes
  all.edges <- data.frame(parent=phy$edge[,1], child=phy$edge[,2])
  
  # extract just the branches that lead to tips
  tip.edges <- subset(all.edges, all.edges$child <= Ntip(phy))
  tip.edges$gCF <- 100 # give these branches arbitrary values of 100
  
  # extract just the internal branches
  int.edges <- subset(all.edges, all.edges$child > Ntip(phy))
  
  # get the branch-appropriate concordance factors
  int.edges <- left_join(int.edges, cf.data)
  int.edges <- int.edges[,c("parent", "child", "gCF")] # subset to just the data we care about
  
  # combine the tip and internal branch information, sort it by the original order!
  combo.edges <- rbind(tip.edges, int.edges)
  combo.edges <- combo.edges[order(match(combo.edges$child, all.edges$child)),]     
  
  # create a color ramp of your choice
  clz <- colorRampPalette(brewer.pal(6, col.palette)); point.colors <- clz(100)
  edge.cols <- point.colors[combo.edges$gCF] # match the CF values to colors
  combo.edges$gCF_col <- edge.cols # add it to the data frame
  
  # make the branches we don't have CFs for black or colored
  if (!is.null(terminal.color)) {
    combo.edges[which(combo.edges$parent == min(combo.edges$parent) & combo.edges$child > Ntip(phy)),]$gCF_col <- terminal.color
    #combo.edges[which(combo.edges$parent == min(combo.edges$parent)+1 & combo.edges$child > Ntip(phy)),]$gCF_col <- "black"
    combo.edges[which(combo.edges$child <= Ntip(phy)),]$gCF_col <- terminal.color
  }

  # if you want to change the line types (it doesn't look great so I'm ignoring this)
  #combo.edges$gCF_lty <- "solid"
  #combo.edges[which(combo.edges$child <= Ntip(phy)),]$gCF_lty <- "dotted"
  
  # provide a check to make sure the values all add up
  if (value.check==TRUE) {tip.list <- combo.edges[which(combo.edges$child <= Ntip(phy)),]
  node.list <- table(tip.list$parent); node.list <- which(node.list > 1)
  target.node <- as.numeric(names(node.list)[sample(1:length(node.list),1)])
  chosen.tips <- getDescendants(phy, target.node)
  #t1 <- phy$tip.label[chosen.tips[1]]; t2 <- phy$tip.label[chosen.tips[2]]; tips <- c(t1,t2)
  true.gCF <- cf.data[which(cf.data$child==target.node),"gCF"]
  target.gCF <- combo.edges[which(combo.edges$child==target.node) ,"gCF"]
  #print(paste("node", target.node, ", the parent of tips", t1, "&", t2, "should have a gCF value of", target.gCF))
  if(!true.gCF==target.gCF){print("something funny about your tree, this plot ain't right")}
  else if(true.gCF==target.gCF){print("all good in the hood")}
  }
  
  # finally, plot it all.
  plot.phylo(phy, edge.color=combo.edges$gCF_col, edge.width=3); # can add 'edge.lty=combo.edges$gCF_lty' 

  if(legend==T){color.bar(colorRampPalette(brewer.pal(6, col.palette))(100))}
  
  layout(mat=matrix(1))
  
  blengths <- data.frame(BranchLength = phy$edge.length, child = phy$edge[,2])
  combo.edges <- left_join(combo.edges, blengths)
  
  return(combo.edges)
}
# e.g.:
# color.by.CF(input.tree, "~/Desktop/cf.stat", terminal.color=NULL, col.palette="RdYlBu", legend=T)

# this function is just to plot the scale bar, which is dumb, but it works, so hey.
color.bar <- function(lut, min=0, max=100, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}