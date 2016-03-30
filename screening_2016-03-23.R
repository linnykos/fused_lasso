rm(list=ls())
source("source_header.R")

registerDoMC(cores = 18)

sigma = 1
n.length = 10
n.vec = round(10^seq(2, 4, length.out = n.length))
trials = 50
jump.mean =     c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)

setup = list(sigma = sigma, n.vec = n.vec, jump.mean = jump.mean, 
             jump.location = jump.location)

level.mat = matrix(0, ncol = trials, nrow = n.length)
haus.haar.mat = matrix(0, ncol = trials, nrow = n.length)
jumps.haar.mat = matrix(0, ncol = trials, nrow = n.length)
haus.filter.mat = matrix(0, ncol = trials, nrow = n.length)
jumps.filter.mat = matrix(0, ncol = trials, nrow = n.length)
haus.filter.boot.mat = matrix(0, ncol = trials, nrow = n.length)
jumps.filter.boot.mat = matrix(0, ncol = trials, nrow = n.length)

simulation_suite <- function(trial){
  set.seed(i*trial*10)
    
  y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)
    
  res = fusedlasso1d(y)
  cv = cv.trendfilter(res, verbose = FALSE)
    
  fit = coef(res, lambda=cv$lambda.1se)$beta
 
  filter.bandwidth = ceiling(0.25*(log(length(y)))^2)
  jumps.filter.idx = apply.filter(fit, filter.bandwidth, 0.5, return.type = "location")
  haus.filter = compute.hausdorff(true.jumps, jumps.filter.idx)
  jumps.filter = length(jumps.filter.idx)

  jumps.haar.idx = apply.filter(y, filter.bandwidth, 0.5, return.type = "location")
  haus.haar = compute.hausdorff(true.jumps, jumps.haar.idx)
  jumps.haar = length(jumps.haar.idx)

  level.fit = bootstrap.threshold(y, fit, filter.bandwidth, cv$lambda.1se,
   trials = 300, quant = 1)
  jumps.filter.boot.idx = apply.filter(fit, filter.bandwidth, level.fit, return.type = "location")
  haus.filter.boot = compute.hausdorff(true.jumps, jumps.filter.boot.idx)
  jumps.filter.boot = length(jumps.filter.boot.idx)

  c(haus.filter, jumps.filter, haus.haar, jumps.haar, level.fit, haus.filter.boot, jumps.filter.boot)
}

for(i in 1:n.length){
  truth = form.truth(jump.mean, jump.location, n.vec[i])
  true.jumps = enumerate.jumps(truth)
 
  res.tmp = foreach(trial = 1:trials) %dopar% simulation_suite(trial)
  res.tmp = do.call(cbind, res.tmp)  

  haus.filter.mat[i,] = res.tmp[1,]
  jumps.filter.mat[i,] = res.tmp[2,]
  haus.haar.mat[i,] = res.tmp[3,]
  jumps.haar.mat[i,] = res.tmp[4,]
  level.mat[i,] = res.tmp[5,]
  haus.filter.boot.mat[i,] = res.tmp[6,]
  jumps.filter.boot.mat[i,] = res.tmp[7,]

  res = list(haus.filter.mat = haus.filter.mat,
   jumps.filter.mat = jumps.filter.mat, haus.haar.mat = haus.haar.mat,
   jumps.haar.mat = jumps.haar.mat, 
   level.mat = level.mat,
   haus.filter.boot.mat = haus.filter.boot.mat,
   jumps.filter.boot.mat = jumps.filter.boot.mat,
   setup = setup)
  save(res, file=paste0("~/DUMP/screening_experiment2_", DATE, ".RData"))
 
  cat('*')
}


