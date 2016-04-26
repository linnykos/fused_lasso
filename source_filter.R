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

bootstrap.threshold <- function(y, fit, filter.bandwidth, lambda, trials = 50, quant = 0.95, verbose = F){
  assert_that(length(y) == length(fit))
  
  n = length(y)
  residuals = y - fit
  
  level.vec = numeric(trials)
  
  for(trial in 1:trials){
    residuals.shuff = residuals[sample(n)]
    
    flasso = fusedlasso1d(residuals.shuff)
    residuals.fit = coef(flasso, lambda = lambda)
    filter.res = apply.filter(residuals.fit$beta, bandwidth = filter.bandwidth, return.type = "filter")

    level.vec[trial] = max(abs(filter.res))

    if(verbose & trial %% floor(trials/10) == 0) cat('*')
  }

  quantile(level.vec, prob = quant)
}

.shuffle.values <- function(y, max.dist = 10){
  d2 = max(floor(max.dist/2), 2)
  
  y2 = y
  n = length(y)

  #first round of permutation
  init.val = sample(d2)[1]
  tmp = init.val  

  y2[1:tmp] = y2[1:tmp][sample(tmp)]
  while(tmp+d2 < n){
    y2[(tmp+1):(tmp+d2)] = y2[(tmp+1):(tmp+d2)][sample(d2)]
    tmp = tmp+d2
  }
  y2[tmp:n] = y2[tmp:n][sample(n-tmp+1)]

  #second round of permutation
  tmp = init.val+d2

  y2[1:tmp] = y2[1:tmp][sample(tmp)]
  while(tmp+d2 < n){
    y2[(tmp+1):(tmp+d2)] = y2[(tmp+1):(tmp+d2)][sample(d2)]
    tmp = tmp+d2
  }
  y2[tmp:n] = y2[tmp:n][sample(n-tmp+1)]

  y2  
}

local.permutation.threshold <- function(y, fit, filter.bandwidth, lambda, trials = 50, 
 quant = 0.9, verbose = F, tol = 1e-3, return.perm = F){
  assert_that(length(y) == length(fit))

  n = length(y)
  filter.org = apply.filter(fit, bandwidth = filter.bandwidth, return.type = "filter")

  custom.func <- function(trial){
    y2 = .shuffle.values(y, filter.bandwidth)
    flasso2 = fusedlasso1d(y2)
    fit2 = coef(flasso2, lambda = lambda)$beta

    apply.filter(fit2, bandwidth = filter.bandwidth, return.type = "filter")
  }

  filter.list = foreach(trial = 1:trials) %dopar% custom.func(trial)

  #now figure out the filter threshold
  #let's do a binary search-type thing

  ##first find the largest filter value ever
  max.filter = max(filter.org)-1e-3
  min.val = max.filter

  high.val = max.filter
  mid.val = max.filter/2
  low.val = 0

  while(high.val - mid.val > tol){
    #first find the changepoints in filter.org
    idx = intersect(which(filter.org >= mid.val), which(abs(diff(fit)) > 1e-5))

    #next, for each value in idx, make sure there is value filter.bandwidth away that is also
    # above the filter at least quant% of the time based on filter.list
    if(length(idx) > 0){ 
      bool = T
    
      for(i in 1:length(idx)){
        idx.window = c(max(1, idx[i]-filter.bandwidth):min(n, idx[i]+filter.bandwidth))
        bool.vec = sapply(filter.list, function(x){
          if(any(x[idx.window] >= mid.val)) return(TRUE) else return(FALSE)
        })
      
        if(sum(bool.vec) <= quant*length(bool.vec)) {bool = F; break()}
      }
    } else bool = F    

    #if all idx passes, that means the filter can be lower
    # else, filter needs to be higher
    if(bool) {
      high.val = mid.val    
    } else {
      low.val = mid.val
    }
    mid.val = (low.val+high.val)/2

    if(verbose) print(mid.val)
  }
 
  if(return.perm) list(val = mid.val, filter.list = filter.list) else mid.val
}
