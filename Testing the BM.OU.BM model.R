library(phytools)
library(geiger)
library(mvtnorm)

source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/New.Models.adapted.from.Slater.2013.R"); ## now source the function from your WD


####
split.vcv <- function(phy, time) {
  
  mat <- vcv(phy)
  mat1 <- mat
  mat2 <- mat
  n <- nrow(mat)
  root <- max(mat)
  shift.from.root <- root - time
  
  make.mat.1 <- function(x, shift.from.root) {
    if(x == 0) {
      return(0)
    }
    if(x<shift.from.root) {
      return(x)
    } 
    if(x > shift.from.root){
      return(shift.from.root)
    }
  }
  
  make.mat.2 <- function(x, time, shift.from.root) {
    if(x == 0) {
      return(0)
    }
    
    if(x < shift.from.root) {
      return(0)
    } else{
      return(x-shift.from.root)
    }
  }
  
  
  mat1 <- matrix(sapply(mat1, make.mat.1, shift.from.root =  shift.from.root), nrow = n, ncol = n, byrow=T)
  
  diag1 <- diag(mat);
  diag.foo <- function(x) {
    if(x<time) {
      return(x)
    } 
    if(x > time){
      
    } 
  }
  mat2 <- matrix(sapply(mat2, make.mat.2, time =  time, shift.from.root = shift.from.root), nrow = n, ncol = n, byrow=T)
  
  rownames(mat1) <- rownames(mat2) <- rownames(mat)
  colnames(mat1) <- colnames(mat2) <- colnames(mat)
  
  return(list(mat1 = mat1, mat2 = mat2))
  
}

####
split.3.vcv <- function(phy, time1, time2) {
  
  mat <- vcv(phy)
  mat1 <- mat
  mat2 <- mat
  mat3 <- mat
  n <- nrow(mat)
  root <- max(mat)
  shift2.from.root <- root - time2
  shift1.from.root <- root - time1
  
  make.mat.1 <- function(x, shift1.from.root) {
    if(x == 0) { #this makes all 0 = 0
      return(0)
    }
    if(x<shift1.from.root) { #this keeps all values that are below s1.f.r (older than the shift1)
      return(x)
    } 
    if(x > shift1.from.root){ #this makes all values above s1.f.r = s1.f.r (equal to the shift1 age)
      return(shift1.from.root)
    }
  }
  
  make.mat.2 <- function(x, time1, time2, shift1.from.root, shift2.from.root) {
    if(x == 0) { #this keeps all 0 = 0
      return(0)
    }
    if((shift2.from.root > x) & (x > shift1.from.root)) { #this makes all values that are above s2.f.r = 0 (makes younger than shift2 = 0)
      return(x)
    }else{
      return(0)
    }
  }
  
  make.mat.3 <- function(x, time2, shift2.from.root) {
    if(x == 0) {
      return(0)
    }
    
    if(x < shift2.from.root) {
      return(0)
    } 
    else{
      return(x-shift2.from.root)
    }
  }
  
  mat1 <- matrix(sapply(mat1, make.mat.1, shift1.from.root =  shift1.from.root), nrow = n, ncol = n, byrow=T)
  
  diag1 <- diag(mat);
  diag.foo <- function(x) {
    if(x<time) {
      return(x)
    } 
    if(x > time){
      
    } 
  }
  mat2 <- matrix(sapply(mat2, make.mat.2, time1 =  time1, time2 = time2, shift1.from.root = shift1.from.root, shift2.from.root = shift2.from.root), nrow = n, ncol = n, byrow=T)
  mat3 <- matrix(sapply(mat3, make.mat.3, time2 =  time2, shift2.from.root = shift2.from.root), nrow = n, ncol = n, byrow=T)
  
  rownames(mat1) <- rownames(mat2) <- rownames(mat3) <- rownames(mat)
  colnames(mat1) <- colnames(mat2) <- rownames(mat3) <- colnames(mat)
  
  return(list(mat1 = mat1, mat2 = mat2, mat3 = mat3))
  
}


#######################################################
# Read in all the data you'll use
######################################################
trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/PB.Meliphagides.100.trees")

data  <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Meliphagides.logMASS.csv", 
                       row.names = 1, header=F) #read in data file in GEIGER format
#data.OUwie <- read.csv("BT.Meliphagides.logMASS.csv", header=F) #read in data file in OUwie format
#data.fitenv<- data.OUwie[,2]; names(data.fitenv) <- data.OUwie[,1] #read in data file in RPANDA format


#logBL<-setNames(data[,14], rownames(data)) #choose your data column (here: logBL) and apply rownames
name.check(trees[[3]], data); #check to make sure the tips match the data labels

test <- vcv(tree)
testvcv <- split.vcv(tree, 10)
testvcv3 <- split.3.vcv(tree, 15,10)
recon <- testvcv3[[1]] + testvcv3[[2]] + testvcv3[[3]]
testvcv3[[1]]
View(testvcv3[[1]])

fitContinuous_paleo(tree, data, model="BM1.OU.BM2", shift.time=5, shift.time2=3)
fitContinuous_paleo(tree, data, model="BM1.OU.BM1", shift.time=7, shift.time2=4)

fitContinuous_paleo(tree, data, model="TRC", shift.time=5)
fitContinuous_paleo(tree, data, model="SRC", shift.time=6)
fitContinuous_paleo(tree, data, model="altTRC", shift.time=5)

fitContinuous(tree, data, model="BM")
fitContinuous_paleo(tree, data, model="BM")

test <- fitContinuous(tree, data, model="OU")
testbnd <- fitContinuous(tree, data, model="OU", bounds=list(alpha=c(0, 50), sigsq=c(0, 10)))
test$bnd



