rm(list=ls())
source("source_header.R")

sigma = 1
n.length = 10
n.vec = round(10^seq(2,4,length.out = n.length))
trials = 50
jump.mean.org =     c(0, 2,  4, 1, 4)
jump.location.org = c(0, .2, .4, .6, .8)

dist.mat = matrix(0, ncol = trials, nrow = n.length)
lambda.vec = matrix(0, ncol = trials, nrow = n.length)

for(i in 1:n.length){
  truth = form.truth(jump.mean.org, jump.location.org, n.vec[i])
  true.jumps = enumerate.jumps(truth)
  
  for(trial in 1:trials){
    set.seed(i*trial*10)
    
    y = generate.problem(n.vec[i], jump.mean.org, jump.location.org,  sigma)
    
    res = fusedlasso1d(y)
    cv = cv.trendfilter(res, verbose = FALSE)
    
    fit = coef(res, lambda=cv$lambda.1se)$beta
    
    dist.mat[i,trial] = compute.hausdorff(true.jumps, enumerate.jumps(fit), 
                                          one.sided = TRUE)
    lambda[i,trial] = cv$lambda.1se
  }
  
  cat('*')
}

setup = list(sigma = sigma, n.vec = n.vec, jump.mean = jump.mean.org, 
             jump.location = jump.location.org)
res = list(dist.mat = dist.mat, lambda.mat = lambda.mat, setup = setup)
save(res, file=paste0("~/DUMP/screening_experiment_", DATE, ".RData"))