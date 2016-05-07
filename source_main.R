##Standard to abide by: a jump at idx "2" means x[2] != x[3]


compute.dual <- function(y, res){
  tmp = res-y
  cumsum(tmp)
}

generate.problem <- function(n, jump.mean, jump.location,  sigma = 0, random.seq = NA){
  #convert jump.location into integer indices
  tmp = seq(0,1,length.out=n)
  jump.location = sapply(jump.location,function(x){min(which(tmp>=x))})
  jump.location[1] = 1
  
  y = rep(0,n)
  for(i in 1:(length(jump.location)-1)){
    tmp = jump.location[i]:(jump.location[i+1]-1)
    y[tmp] = rnorm(length(tmp),mean = jump.mean[i], sd = sigma)
  }
  tmp = jump.location[length(jump.location)]:n
  y[tmp] = rnorm(length(tmp),mean = jump.mean[length(jump.mean)], sd = sigma)
  
  if(!missing(random.seq)){
    y = y+random.seq
  }
  
  y
}


dualplot_suite <- function(y, jump.mean, jump.location, lambda = NA, 
                           num.est.jumps = NA, plot.res = TRUE,
                           fused.lasso.obj = NA, return.fusedlasso = FALSE){
  assert_that(!is.na(lambda) || !is.na(num.est.jumps))
  
  if(is.na(fused.lasso.obj[1])) {
    fres = fusedlasso1d(y)
  } else {
    fres = fused.lasso.obj
  }
  
  if(!is.na(lambda)) {
    tmp = coef(fres,lambda=lambda)
    res = tmp$beta
  } else {
    tmp = coef(fres,df=num.est.jumps+1)
    res = tmp$beta
    lambda = tmp$lambda
  }

  ret = summarize(res, jump.mean, jump.location, lambda)
  
  if(plot.res) plotfused(jump.mean,jump.location,y,ret[["coef"]],lambda,ret[["jumps"]],ret[["mse"]])
  
  ret
}

summarize <- function(fit, jump.mean, jump.location, lambda){
  mse = compute.mse(fit, jump.mean, jump.location)
  num.est.jumps = count.jumps(fit)
  
  ret = list()
  class(ret) = "fusedLassoSummary"
  ret[["coef"]] = fit
  ret[["lambda"]] = lambda
  ret[["jumps"]] = num.est.jumps
  ret[["mse"]] = mse
  
  return(ret)
}


compute.mse <- function(res, jump.mean=NA, jump.location=NA, true.seq=NA){
  bool1 = !is.na(jump.mean[1]) & !is.na(jump.location[1])
  bool2 = !is.na(true.seq[1])
  assert_that(bool1 | bool2)
  
  n = length(res)
  
  if(bool1) true.seq = form.truth(jump.mean, jump.location, n)
  
  sum((true.seq-res)^2)/n
}

count.jumps <- function(fit, tol=1e-3){
  length(enumerate.jumps(fit, tol))
}

form.truth <- function(jump.mean, jump.location, n){
  assert_that(min(jump.location)>=0)
  assert_that(max(jump.location)<=1)
  
  tmp = seq(0,1,length.out=n)
  jump.location2 = sapply(jump.location,function(x){max(min(which(tmp>=x)),1)-1})
  jump.location2[1] = 1
  jump.location2 = sort(jump.location2)
  jump.location3 = c(jump.location2,n)
  jump.location3[1] = 0
  
  rep(jump.mean, times=diff(jump.location3))
}

#outputs the first index in a new segment
enumerate.jumps <- function(fit, tol = 1e-4, include.endpoints = FALSE){
  dif = abs(diff(fit))
  idx = which(dif > tol) + 1
  
  if(include.endpoints) c(1, idx, length(fit)+1) else idx
}

#if one.sided = TRUE, then only measure the distance from set1 to set2
## (i.e., how well set2 covers set1)
compute.hausdorff <- function(set1, set2, one.sided = FALSE){

  #handle corner cases
  if(length(set1) == 0 | length(set2) == 0) return(NA)
  if(length(set2) == 1) set2 = c(set2, set2)
  if(length(set1) == 1) set1 = c(set1, set1)

  dist.mat = sapply(set1, function(i){abs(i-set2)})
  if(class(dist.mat) != "matrix") dist.mat = as.matrix(dist.mat, nrow = 1)
  dist.vecx = apply(dist.mat, 2, min)
  
  if(!one.sided) dist.vecy = apply(dist.mat, 1, min) else dist.vecy = 0
  
  max(dist.vecx, dist.vecy)
}

form.diff.matrix <- function(n){
  res = matrix(0,n-1,n)
  for(i in 1:(n-1)){
    res[i,c(i,i+1)] = c(-1,1)
  }
  
  res
}
