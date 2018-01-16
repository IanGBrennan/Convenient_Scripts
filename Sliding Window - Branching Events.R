library(evobiR)
library(phytools)
library(laser)
rainbowfish<-read.tree("oz.rainbowfish.tre")
plot(rainbowfish)
rfish.btimes<-branching.times(rainbowfish)

yuleWindow(rfish.btimes, 70, 65)

sort(rfish.btimes, decreasing=TRUE)

SlidingWindow(yuleWindow, rfish.btimes, 5, 1)

SlidingWindow
yuleWindow


x <- rev(sort(rfish.btimes))
res <- as.list(yuleint2(rfish.btimes, 49, 47))
st1<-49
st2<-47



data.rf<-rfish.btimes
# the embedded function, modified from yuleWindow in LASER
crackedWindow<-function (data.rf, st.x, st.y, step) 
{
  nvec <- 2:(length(data.rf) + 1)
  nv <- data.rf[(data.rf < st.x) & (data.rf >= st.y)]
  lo <- max(nvec[data.rf >= st.x])
  up <- max(nvec[data.rf >= st.y])
  spots <- seq(from=1, to=max(data.rf), by=1)
  result<-vector(length=length(spots))
  if (st.x <= data[1]) {
    nv <- c(st.x, nv) - st.y
  }
  else {
    nv <- nv - st.y
  }
  b.int<-(up-lo)/lo
  for (i in 1:length(spots)) {
    result[i]<-b.int
  }
}


crackedWindow<-function (data, st.x, st.y) 
{
  nvec <- 2:(length(data) + 1)
  nv <- data[(data < st.x) & (data >= st.y)]
  lo <- max(nvec[data >= st.x])
  up <- max(nvec[data >= st.y])
  if (st.x <= data[1]) {
    nv <- c(st.x, nv) - st.y
  }
  else {
    nv <- nv - st.y
  }
  b.int<-(up-lo)/lo
  st.x
}

est<-crackedWindow(rfish.btimes, 5, 0)

# the parent function, modified from SlidingWindow in evobiR 
slidingDoor<-function (FUN, data, window, step) 
{
  total <- max(data)
  spots <- seq(from = 1, to = (total - window), by = step)
  result <- vector(length = length(spots))
  for (i in 1:length(spots)) {
    result[i] <- match.fun(FUN)(data[spots[i]:(spots[i] + window-1)])
  }
  return(result)
}
st.1<-max(rfish.btimes)
st.1<-45
first<-crackedWindow(rfish.btimes, 5, 0)
#repeat for whole of tree
