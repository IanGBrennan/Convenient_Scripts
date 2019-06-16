# Haploid population size
N = 10
# Number of sampled gene copies
k = 6
# Upper limit for the plot
#tlim = 2*N*((k-1)/N)
tlim = 30

# Delay for the animation (in seconds)
delay <- 0.05

freeze <- function(){ Sys.sleep(delay)}
# Sampling time
t= 0

# Plot the haploid individuals of the whole population
pop <- seq(from = 1, to = N, by = 1)
plot(pop, rep(t,N), ylim = c(0,tlim), xlab = "individuals", ylab = "time")

# Simulate in black the sampled gene copies
sample <- sample(x = pop, size = k)
points(sample, rep(t,k), pch=16)

coal.times <- NULL
tmrca <- NULL
# Loop to coalesce until the MRCA (or to stop when outside the plot)
while(k > 1 && t < tlim){
  t = t + 1
  # Plot the previous population of individuals
  points(pop, rep(t,N))
  freeze()
  # Just to store the parents that have been picked by the lineages
  remember <- c()
  # Iterating over the lineages
  for(child in sort(sample)){ #child = sort(sample)[1]
    # Sample one parent uniformely at random
    parent = sample(x = pop, size = 1)
    # Plot the ancestry relationship between the child and the parent
    #segments(x0 = child, y0 = t-1, x1 = parent, y1 = t)
    freeze()
    # Is it a coalescence event ?
    number_of_coalescence = sum(remember == parent)
    if(number_of_coalescence == 0){
      # Nope :(
      points(parent, t, pch=16, col="black")
    }else if (number_of_coalescence == 1){
      # Yes ! :D
      points(parent, t, pch=16, col="blue", cex=2)
      print(paste("coalescent event at time:", t))
      coal.times <- append(coal.times, t)
      k = k-1
    }else if (number_of_coalescence >= 2){
      # OMG more than two nodes have coalesced!!!!
      points(parent, t, pch=16, col="red", cex=2)
      print(paste("double coalescent event at time", t))
      coal.times <- append(coal.times, t)
      k = k-1
    }
    freeze()
    # Add the parent to the list of already picked parents
    remember<-c(remember,parent)
  }
  # The remaining lineages in the parental generation form the sample of the next iteration.
  sample<-unique(remember)
  if(length(sample)==1){
    tmrca = t
    points(parent, t, pch=16, col="green", cex=2)
    print(paste("TMRCA reached at time", t))}
}
#plot(density(coal.times))
