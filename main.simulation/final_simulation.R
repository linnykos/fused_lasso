#load sourcese
rm(list=ls())
setwd("~/ryan/fused.git/")
source("source_header.R")

#set up parameters
sigma = 2
n.length = 10
n.vec = round(10^seq(2,4,length.out=n.length))
trials = 50
jump.mean = c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)
haus.quant = 0.95

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
  threshold.oracle = oracle.threshold(fit, filter.bandwidth, truth)

  jumps.oracle = apply.filter(fit, filter.bandwidth, threshold.oracle, 
   return.type = "location")
  left.oracle = compute.hausdorff(true.jumps, jumps.oracle, one.sided = T)
  right.oracle = compute.hausdorff(jumps.oracle, true.jumps, one.sided = T)

  #run the filter
  threshold.adapt = staircase.threshold(y, fit, filter.bandwidth, cv$lambda.min,
   controls = list(type = "original", quant = seq(0, 1, length.out = 101)))

  idx = which(names(threshold.adapt) == paste0(100*haus.quant, "%"))

  jumps.adapt = apply.filter(fit, filter.bandwidth, threshold.adapt[idx],
   return.type = "location")
  left.adapt = compute.hausdorff(true.jumps, jumps.adapt, one.sided = T)
  right.adapt = compute.hausdorff(jumps.adapt, true.jumps, one.sided = T)

  list(lambda = cv$lambda.min, mse = mse,
   left.org = left.org, right.org = right.org,
   left.oracle = left.oracle, right.oracle = right.oracle,
   left.adapt = left.adapt, right.adapt = right.adapt,
   threshold.oracle = threshold.oracle, threshold.adapt = threshold.adapt)
}

#set up the matrices where we're going to store our results
res.list = lapply(1:9, function(x){
  mat = matrix(0, ncol = n.length, nrow = trials)
  colnames(mat) = n.vec
  mat
})
res.list[[10]] = lapply(1:n.length, function(x){
  mat = matrix(0, ncol = 101, nrow = trials)
  colnames(mat) = seq(0, 1, length.out = 101)
  mat
})

names(res.list) = c("lambda", "mse", "left.org", "right.org", "left.oracle",
  "right.oracle", "left.adapt", "right.adapt", "threshold.oracle", 
  "threshold.adapt")
names(res.list[[10]]) = n.vec

#run the simulations
for(i in 1:length(n.vec)){
  truth = form.truth(jump.mean, jump.location, n.vec[i])
  
  for(trial in 1:trials){
    set.seed(i*trial*10)
    y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)

    res = run.tests(y, truth)

    #unpack the results
    for(j in 1:(length(res.list)-1)){
      res.list[[j]][trial,i] = res[[j]]
    }
    res.list[[length(res.list)]][[i]][trial,] = res[[length(res.list)]]

    save.image(file = paste0("results/final-", Sys.Date(), ".RData"))
  }
}


