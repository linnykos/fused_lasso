#more of just simple experiments than tests
rm(list=ls())
setwd("..")
source("source_header.R")

sigma = 3.5
n.length = 1
n.vec = 1000
trials = 25
jump.mean =     c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)

registerDoMC(cores = 20)

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
   
     print(paste0("Starting on ",  names(filter.list)[j]))

     set.seed(i*trial*j*10)
 
     filter.list[[j]][[i]][trial,] = staircase.threshold(y, fit, 
      filter.bandwidth, cv$lambda.1se, controls = list(type = 
      names(filter.list)[j], quant = seq(0, 1, length.out = 
      ncol(filter.list[[j]][[i]]))))
  
     res = list(filter.list = filter.list, setup = setup)
     save(res, file = paste0("~/ryan/fused.git/results/filterExperiment_", 
      DATE, ".RData"))
    }

    print(paste0("Trial ", trial, " for i=", i, " complete!"))
  }
}

