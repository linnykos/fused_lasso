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

  jumps.org = enumerate.jumps(fit)

  jumps.oracle.list = vector("list", length(oracle.seq))
  jumps.adapt.list = vector("list", length(oracle.seq))

  for(k in 1:length(oracle.seq)){
    jumps.oracle.list[[k]] = apply.filter(fit, filter.bandwidth,
     oracle.seq[k], return.type = "location")

    jumps.adapt.list[[k]] = apply.filter(fit, filter.bandwidth,
     filter.mat[trial,k], return.type = "location")
  }

  list(jumps.org = jumps.org, jumps.oracle.list = jumps.oracle.list, 
   jumps.adapt.list = jumps.adapt.list)
}

res.list = vector("list", 3)
names(res.list) = c("jumps.org", "jumps.oracle.list", "jumps.adapt.list")
res.list[[1]] = vector("list", trials)
for(k in 2:3){
  res.list[[k]] = vector("list", length(oracle.seq))
  for(j in 1:length(oracle.seq)) res.list[[k]][[j]] = vector("list", trials)
}

#run the simulations
truth = form.truth(jump.mean, jump.location, n.vec[i])

for(trial in 1:trials){
  set.seed(i*trial*10)
  y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)

  res = run.test(y, truth, trial)

  #unpack the results
  res.list[[1]][[trial]] = res$jumps.org
  for(k in 1:length(oracle.seq)){
    res.list[[2]][[k]][[trial]] = res$jumps.oracle.list[[k]]
    res.list[[3]][[k]][[trial]] = res$jumps.adapt.list[[k]]
  }

  save.image(file = paste0("~/ryan/fused.git/results/final-ROC2-", Sys.Date(), ".RData"))
 
}


