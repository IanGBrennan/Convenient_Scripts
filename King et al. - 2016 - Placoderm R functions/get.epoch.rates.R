# Can this script be adjusted to pull the rates of biome transition out?
# then we can compare the rates of biome movement through time periods,
# and investigate a speed up in the Miocene

get.epoch.rates <- function(trees, epochs, burnin){
	epochs-epochs[length(epochs)] -> binages
	data.frame(binmin=binages[2:length(epochs)], binmax=binages[1:length(epochs)-1]) -> epochbounds
	trees[((round(length(trees)*burnin))+1):length(trees)] -> treesb
	matrix(0, length(treesb), length(epochs)-1) -> meansample
	for(i in 1:(length(treesb))){
		node.age(treesb[[i]]) -> nodeages
		sapply(nodeages$ages, function(x) max(nodeages$ages) - x) -> minheights
		treesb[[i]]$edge.length + minheights -> maxheights
		data.frame(minheight=minheights, maxheight=maxheights, rate=1:length(treesb[[1]]$edge.length)) -> branchrates
		for(j in 1:(length(treesb[[1]]$edge.length))){
			treesb[[i]]$annotations[[j]]$rate -> branchrates[j,3]
		}
		matrix(0, length(treesb[[1]]$edge.length), length(epochs)-1) -> epochbranchlens
		for(j in 1:((length(epochs))-1)){
			for(k in 1:(length(treesb[[1]]$edge.length))){
				if(branchrates[k,1] >= epochbounds[j,2])
				0 -> epochbranchlens[k,j]
				if(branchrates[k,1] < epochbounds[j,2])
					if(branchrates[k,1] >= epochbounds[j,1])
						if(branchrates[k,2] >= epochbounds[j,2])
						epochbounds[j,2] - branchrates[k,1] -> epochbranchlens[k,j]
				if(branchrates[k,1] < epochbounds[j,1])
					if(branchrates[k,2] >= epochbounds[j,2])
					epochbounds[j,2] - epochbounds[j,1] -> epochbranchlens[k,j]
				if(branchrates[k,2] < epochbounds[j,2])
					if(branchrates[k,1] >= epochbounds[j,1])
					branchrates[k,2] - branchrates[k,1] -> epochbranchlens[k,j]
				if(branchrates[k,2] < epochbounds[j,2])
					if(branchrates[k,2] >= epochbounds[j,1])
						if(branchrates[k,1] < epochbounds[j,1])
						branchrates[k,2] - epochbounds[j,1] -> epochbranchlens[k,j]
				if(branchrates[k,2] < epochbounds[j,1])
					0 -> epochbranchlens[k,j]
			}
		}
		epochbranchlens -> epochratesweight
		for(j in 1:(length(treesb[[1]]$edge.length))){
			epochratesweight[j,1:((length(epochs))-1)]*branchrates[j,3] -> epochratesweight[j,1:((length(epochs))-1)]
		}
		apply(epochratesweight, 2, sum)/apply(epochbranchlens, 2, sum) -> epochmeans
		epochmeans -> meansample[i,]
	}
	return(meansample)
}

trees<-read.nexus("Ag100_3clades_strGTI_LinkedMuGTI_all.trees")
epochs<-c(40, 30, 20, 10, 0)
ratez<-get.epoch.rates(trees, epochs, burnin=0)

