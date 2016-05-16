#load sourcese
rm(list=ls())
setwd("~/ryan/fused.git/")
source("source_header.R")

load("~/ryan/fused.git/results/final-ROC-2016-05-13.RData")

registerDoMC(cores = 20)

oracle.seq = seq(0, 2, length.out = ncol(filter.mat))

#testing function
run.test <- function(y, truth, trial){
  res = fusedlasso1d(y)
  cv = cv.trendfilter(res, verbose = FALSE)

  fit = coef(res, lambda=cv$lambda.min)$beta

  n = length(y)
  filter.bandwidth = ceiling(0.25*(log(n))^2)
  true.jumps = enumerate.jumps(truth)

  oracle.left = numeric(length(oracle.seq))
  oracle.right = oracle.left
  adapt.left = oracle.left
  adapt.right = oracle.left

  for(k in 1:length(oracle.seq)){
    jumps = apply.filter(fit, filter.bandwidth,
     oracle.seq[k], return.type = "location")
    oracle.left[k] = compute.hausdorff(true.jumps,
     jumps, one.sided = T)
    oracle.right[k] = compute.hausdorff(jumps,
     true.jumps, one.sided = T)

    jumps = apply.filter(fit, filter.bandwidth,
     filter.mat[trial,k], return.type = "location")
    adapt.left[k] = compute.hausdorff(true.jumps,
     jumps, one.sided = T)
    adapt.right[k] = compute.hausdorff(jumps,
     true.jumps, one.sided = T)
  }

  list(oracle.left = oracle.left, oracle.right = oracle.right,
   adapt.left = adapt.left, adapt.right = adapt.right)
}

res.list = vector("list", 4)
names(res.list) = c("oracle.left", "oracle.right", "adapt.left", "adapt.right")
for(k in 1:4){
  res.list[[k]] = matrix(0, ncol = length(oracle.seq), nrow = trials)

  if(k == 1 | k == 2) colnames(res.list[[k]]) = oracle.seq else
    colnames(res.list[[k]]) = paste0(colnames(filter.mat), "%")
}

#run the simulations
truth = form.truth(jump.mean, jump.location, n.vec[i])

for(trial in 1:trials){
  set.seed(i*trial*10)
  y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)

  res = run.test(y, truth, trial)

  #unpack the results
  for(k in 1:4) res.list[[k]][trial,] = res[[k]]

  save.image(file = paste0("~/ryan/fused.git/results/final-ROC2-", DATE, ".RData"))
 
}


