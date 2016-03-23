apply.filter <- function(fit, bandwidth, threshold = 0, y = NA, return.type = c("location", "filter", "refit")){
  n = length(fit)
  return.type = return.type[1]
  
  if(return.type == "refit") assert_that(all(!is.na(y)))
  
  jump.loc = enumerate.jumps(fit)

  z = sapply((bandwidth+1):(n-bandwidth+1), function(x){
    mean(fit[x:(x+bandwidth-1)]) - mean(fit[(x-bandwidth):(x-1)])
  })
  
  #pad the z's
  z = c(rep(0, bandwidth), z, rep(0, bandwidth-1))
  if(return.type == "filter") return(z)
  
  z.thresIdx = which(abs(z) >= threshold)
  
  jump.remain = intersect(jump.loc, z.thresIdx)
  if(return.type == "location") return(jump.remain)
  
  jump.remain = c(1, jump.remain, n+1)
  refit = sapply(2:length(jump.remain), function(x){
    mean(y[x-1]:(y[x]-1))
  })
  
  return(refit)
}

bootstrap.threshold <- function(y, fit, filter.bandwidth, lambda, trials = 500, quant = 0.95){
  assert_that(length(y) == length(fit))
  
  n = length(y)
  residuals = y - fit
  
  level.vec = numeric(trials)
  
  for(trial in 1:trials){
    residuals.shuff = residuals[sample(n)]
    
    flasso = fusedlasso1d(residuals.shuff)
    residuals.fit = coef(flasso, lambda = lambda)
    filter.res = apply.filter(flasso, bandwidth = filter.bandwidth, return.type = "filter")

    level.vec[trial] = max(abs(filter.res))
  }

  quantile(level.vec, prob = quant)
}

