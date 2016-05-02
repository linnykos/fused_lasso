#more of just simple experiments than tests
rm(list=ls())
setwd("..")
source("source_header.R")

sigma = 1
n.length = 10
n.vec = round(10^seq(2, 4, length.out = n.length))
trials = 50
jump.mean =     c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)

registerDoMC(cores = 10)

setup = list(sigma = sigma, n.vec = n.vec, jump.mean = jump.mean, 
  jump.location = jump.location)

filter.list = vector("list", 3)
for(j in 1:3){
  filter.list[[j]] = vector("list", n.length)

  for(i in 1:n.length){
    filter.list[[j]][[i]] = matrix(0, ncol = 101, nrow = trials)
  }

  names(filter.list[[j]]) = n.vec
}
names(filter.list) = c("original", "rearrange", "flip")

for(i in 1:length(n.vec)){
  truth = form.truth(jump.mean, jump.location, n.vec[i])
  true.jumps = enumerate.jumps(truth)

  filter.bandwidth = ceiling(0.25*(log(length(truth)))^2)

  for(trial in 1:trials){
    set.seed(i*trial*10)

    y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)
    
    res = fusedlasso1d(y)
    cv = cv.trendfilter(res, verbose = FALSE)
    
    fit = coef(res, lambda=cv$lambda.1se)$beta
 
    for(j in 1:length(filter.list)){

     filter.list[[j]][[i]][trial,] = staircase.threshold(y, fit, 
      filter.bandwidth, cv$lambda.1se, controls = list(type = 
      names(filter.list)[j], quant = seq(0, 1, length.out = 
      ncol(filter.list[[j]][[i]]))))
  
     res = list(filter.list = filter.list, setup = setup)
     save(res, file = paste0("~/ryan/fused.git/results/filterExperiment_", 
      DATE, ".RData"))
    }
  }
}

