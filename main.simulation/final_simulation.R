#load sourcese
rm(list=ls())
setwd("~/ryan/fused.git/")
source("source_header.R")

#set up parameters
sigma = 3.5
n.length = 10
n.vec = round(10^seq(2,4,length.out=n.length))
trials = 25
jump.mean = c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)
haus.quant = 0.8

registerDoMC(cores = 20)

#set up the function the run all our tests for us
run.tests <- function(y, truth){
  res = fusedlasso1d(y)
  cv = cv.trendfilter(res, verbose = F)

  fit = coef(res, lambda = cv$lambda.min)$beta
  mse = compute.mse(fit, true.seq = truth)

  #reconstruct the parameters
  n = length(y)
  filter.bandwidth = ceiling(0.25*(log(n))^2)
  true.jumps = enumerate.jumps(truth)

  #compute the hausdorff without filter
  jumps.org = enumerate.jumps(fit)
  left.org = compute.hausdorff(true.jumps, jumps.org, one.sided = T)
  right.org = compute.hausdorff(jumps.org, true.jumps, one.sided = T)

  #run the oracle filter
  jumps.oracle = apply.filter(fit, filter.bandwidth, 0.5, 
   return.type = "location")
  left.oracle = compute.hausdorff(true.jumps, jumps.oracle, one.sided = T)
  right.oracle = compute.hausdorff(jumps.oracle, true.jumps, one.sided = T)

  #run the filter
  threshold = staircase.threshold(y, fit, filter.bandwidth, cv$lambda.min,
   controls = list(type = "original", quant = 0.8))
  jumps.adapt = apply.filter(fit, filter.bandwidth, threshold,
   return.type = "location")
  left.adapt = compute.hausdorff(true.jumps, jumps.adapt, one.sided = T)
  right.adapt = compute.hausdorff(jumps.adapt, true.jumps, one.sided = T)

  list(lambda = cv$lambda.min, mse = mse,
   left.org = left.org, right.org = right.org,
   left.oracle = left.oracle, right.oracle = right.oracle,
   threshold = threshold, left.adapt = left.adapt, right.adapt = right.adapt)
}

#set up the matrices where we're going to store our results
res.list = lapply(1:9, function(x){
  mat = matrix(0, ncol = n.length, nrow = trials)
  colnames(mat) = n.vec
  mat
})
names(res.list) = c("lambda", "mse", "left.org", "right.org", "left.oracle",
  "right.oracle", "threshold", "left.adapt", "right.adapt")

#run the simulations
for(i in 1:length(n.vec)){
  truth = form.truth(jump.mean, jump.location, n.vec[i])
  
  for(trial in 1:trials){
    set.seed(i*trial*10)
    y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)

    res = run.tests(y, truth)

    #unpack the results
    for(j in 1:length(res.list)){
      res.list[[j]][trial,i] = res[[j]]
    }

    save.image(file = paste0("~/ryan/fused.git/results/final-", DATE, ".RData"))
  }
}


