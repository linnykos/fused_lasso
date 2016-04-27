#more of just simple experiments than tests
rm(list=ls())
setwd("..")
source("source_header.R")

sigma = 0.2
n.length = 10
n.vec = round(10^seq(2, 4, length.out = n.length))
trials = 50
jump.mean =     c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)

registerDoMC(cores = 10)

i = 4
truth = form.truth(jump.mean, jump.location, n.vec[i])
true.jumps = enumerate.jumps(truth)

setup = list(sigma = sigma, n.vec = n.vec, jump.mean = jump.mean, 
             jump.location = jump.location, n.idx = i)

filter.mat = matrix(0, ncol = trials, nrow = 4)
rownames(filter.mat) = c("Bootstrap.Residual", "Local.Permutation", 
 "Bootstrap.Staircase", "Oracle")

haus.mat = filter.mat
jumps.mat = filter.mat

for(trial in 1:trials){
  set.seed(i*trial*10)

  y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)
    
  res = fusedlasso1d(y)
  cv = cv.trendfilter(res, verbose = FALSE)
    
  fit = coef(res, lambda=cv$lambda.1se)$beta
 
  filter.bandwidth = ceiling(0.25*(log(length(y)))^2)

  #first try the standard bootstrap idea
  filter.mat[1,trial] = bootstrap.threshold(y, fit, filter.bandwidth, cv$lambda.1se,
   trials = 100, quant = 0.95)
  filter.mat[2,trial] = local.permutation.threshold(y, fit, filter.bandwidth, cv$lambda.1se,
   trials = 100, quant = 0.95)
  filter.mat[3,trial] = staircase.threshold(y, fit, filter.bandwidth, cv$lambda.1se,
   trials = 100, quant = 0.95)
  filter.mat[4,trial] = 0.5

  #calculate the haus.mat and jumps.mat
  for(type in 1:nrow(filter.mat)){
    jumps.filter.boot.idx = apply.filter(fit, filter.bandwidth, 
     filter.mat[type,trial], return.type = "location")
    haus.mat[type,trial] = compute.hausdorff(true.jumps, jumps.filter.boot.idx)
    jumps.mat[type,trial] = length(jumps.filter.boot.idx)
  }

  res = list(filter.mat = filter.mat, haus.mat = haus.mat, jumps.mat = jumps.mat,
   setup = setup)
  save(res, file=paste0("~/ryan/fused.git/results/filterExperiment_n-", n.vec[i], 
   "-", DATE, ".RData"))
}


