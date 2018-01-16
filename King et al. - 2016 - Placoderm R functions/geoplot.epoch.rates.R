geoplot.epoch.rates <- function(rate.sample, epochs, line.col="black", plot.mean=T, mean.col="red", plot.bounds=T, ...){
	data.frame(binmin=epochs[2:length(epochs)], binmax=epochs[1:((length(epochs))-1)]) -> epochbounds
	apply(data.frame(binmin=epochs[2:length(epochs)], binmax=epochs[1:length(epochs)-1]), 1, mean) -> midages
	geoscalePlot(ages=(midages), data=rate.sample[1,], data.lim=c(min(rate.sample), max(rate.sample)), type="n", label="rate", ...)
	if(plot.bounds == TRUE)
		for(i in 1:length(midages)){
			abline(v=epochbounds[1:length(midages),1][i], lty=3)
		}
	for(i in 1:nrow(rate.sample)){
		lines(midages[1:length(midages)], rate.sample[i,1:length(midages)], col=line.col, lwd=0.05)
	}
	if(plot.mean == TRUE)
		lines(midages[1:length(midages)], apply(rate.sample, 2, mean), col=mean.col, lwd=3)
	
}


