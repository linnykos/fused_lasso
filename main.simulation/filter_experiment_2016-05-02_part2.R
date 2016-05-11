#more of just simple experiments than tests
rm(list=ls())
setwd("~/ryan/fused.git")
source("source_header.R")

sigma = 4
n.length = 1
n.vec = 5000
trials = 50
jump.mean = seq(0, 5, length.out = 6)
jump.location = seq(0, 1, length.out = 7)[-7]

registerDoMC(cores = 16)

setup = list(sigma = sigma, n.vec = n.vec, jump.mean = jump.mean, 
  jump.location = jump.location)

left.list = vector("list", 3)
for(j in 1:3){
  left.list[[j]] = vector("list", n.length)

  for(i in 1:n.length){
    left.list[[j]][[i]] = matrix(0, ncol = 101, nrow = trials)
    colnames(left.list[[j]][[i]]) = paste0(seq(0,1,length.out=101)*100,"%")
  }

  names(left.list[[j]]) = n.vec
}
names(left.list) = c("original", "rearrange", "flip")
right.list = left.list

#creat the oracle lists
oracle.seq = seq(0, 3, length.out = 151)

oracle.left.list = vector("list", n.length)
for(i in 1:n.length){
  oracle.left.list[[i]] = matrix(0, ncol = length(oracle.seq), nrow = trials)
  colnames(oracle.left.list[[i]]) = as.character(oracle.seq)
}
oracle.right.list = oracle.left.list

#load in the thresholds
load("~/ryan/fused.git/results/filterExperiment_2016-05-07.RData")
filter.list = res$filter.list

#run the simulations to see how the hausdorff distance is
#WARNING: HAS BEEN CHANGED
#for(i in 1:length(n.vec)){
for(i in 1:length(n.vec)){
  truth = form.truth(jump.mean, jump.location, n.vec[i])
  true.jumps = enumerate.jumps(truth)

  filter.bandwidth = ceiling(0.25*(log(length(truth)))^2)

  custom.func <- function(trial){
    print(trial)
    set.seed(i*trial*10)

    y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)
    
    res = fusedlasso1d(y)
    cv = cv.trendfilter(res, verbose = FALSE)
    
    fit = coef(res, lambda = cv$lambda.1se)$beta

    left.dist = matrix(0, ncol = 101, nrow = 3)
    right.dist = left.dist

    #evaluate the different filters
    for(j in 1:3){
      for(k in 1:101){
        if(is.na(filter.list[[j]][[i]][trial,k])){
          left.dist[j,k] = NA; right.dist[j,k] = NA
        }

        jumps.filter.boot.idx = apply.filter(fit, filter.bandwidth, 
         filter.list[[j]][[i]][trial,k], return.type = "location")
        left.dist[j,k] = compute.hausdorff(true.jumps, 
         jumps.filter.boot.idx, one.sided = T)
        right.dist[j,k] = compute.hausdorff(jumps.filter.boot.idx,
         true.jumps, one.sided = T)
      }
    }

    #now compute the oracle one
    oracle.left.dist = numeric(151)
    oracle.right.dist = oracle.left.dist

    for(k in 1:151){
      jumps.filter.boot.idx = apply.filter(fit, filter.bandwidth,
       oracle.seq[k], return.type = "location")
      oracle.left.dist[k] = compute.hausdorff(true.jumps,
       jumps.filter.boot.idx, one.sided = T)
      oracle.right.dist[k] = compute.hausdorff(jumps.filter.boot.idx,
       true.jumps, one.sided = T)
    }

    list(right.dist = right.dist, left.dist = left.dist, 
     oracle.left.dist = oracle.left.dist, 
     oracle.right.dist = oracle.right.dist)
  }

  res = foreach(trial = 1:trials) %dopar% custom.func(trial)

  #unpack the foreach results
  for(trial in 1:trials){
    for(j in 1:3){
      left.list[[j]][[i]][trial,] = res[[trial]]$left.dist[j,]
      right.list[[j]][[i]][trial,] = res[[trial]]$right.dist[j,]
    }

    oracle.left.list[[i]][trial,] = res[[trial]]$oracle.left.dist
    oracle.right.list[[i]][trial,] = res[[trial]]$oracle.right.dist
  }

  res = list(left.list = left.list, right.list = right.list, 
   oracle.left.list = oracle.left.list, 
   oracle.right.list = oracle.right.list)
  save(res, file = 
   paste0("~/ryan/fused.git/results/filterExperiment_hausdorff_", 
   Sys.Date(), ".RData"))
}


#create the average-performance ROC curves

