apply.filter <- function(fit, bandwidth, threshold = 0, y = NA, return.type = c("location", "filter", "refit")){
  n = length(fit)
  return.type = return.type[1]
  
  if(return.type == "refit") assert_that(all(!is.na(y)))
  
  jump.loc = enumerate.jumps(fit)

  z = sapply((bandwidth+1):(n-bandwidth+1), function(x){
    abs(mean(fit[x:(x+bandwidth-1)]) - mean(fit[(x-bandwidth):(x-1)]))
  })
  
  #pad the z's
  z = c(rep(0, bandwidth), z, rep(0, bandwidth-1))
  if(return.type == "filter") return(z)
  
  z.thresIdx = which(abs(z) >= threshold)
 
  #WARNING: check this
  jump.loc2 = c(jump.loc, jump.loc-bandwidth, jump.loc+bandwidth-1)
  jump.loc2 = jump.loc2[which(jump.loc2>0)]
  jump.loc2 = jump.loc2[which(jump.loc2<=n)]
  jump.loc2 = sort(unique(jump.loc2))
 
  jump.remain = intersect(jump.loc2, z.thresIdx)
  if(return.type == "location") return(jump.remain)
  
  #WARNING: fix this so its a refit of fused lasso
  jump.remain = c(1, jump.remain, n+1)
  refit = sapply(2:length(jump.remain), function(x){
    mean(y[x-1]:(y[x]-1))
  })
  
  return(refit)
}

#bootstrap the residuals, used fused lasso on the residuals and record filter value to kill the jumps
bootstrap.threshold <- function(y, fit, filter.bandwidth, lambda, trials = 50, quant = 0.95){
  assert_that(length(y) == length(fit))
  
  n = length(y)
  residuals = y - fit
  
  custom.func <- function(trial){
    set.seed(10*trial)
    residuals.shuff = residuals[sample(n)]
    
    flasso = fusedlasso1d(residuals.shuff)
    residuals.fit = coef(flasso, lambda = lambda)
    filter.res = apply.filter(residuals.fit$beta, bandwidth = filter.bandwidth, return.type = "filter")

    max(abs(filter.res))
  }

  level.vec = foreach(trial = 1:trials) %dopar% custom.func(trial)

  quantile(unlist(level.vec), prob = quant)
}
