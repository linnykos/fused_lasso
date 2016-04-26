#rearrange the fit to be a big staircase
staircase.threshold <- function(y, fit, filter.bandwidth, lambda, trials = 50,
  quant = 0.9, tol = 1e-6){

  assert_that(length(y) == length(fit))

  n = length(y)
  residual = y - fit

  staircase = .arrange.staircase(fit, tol)
  sc.jumps = enumerate.jumps(staircase)

  sc.jumpNeigh = sapply(sc.jumps, function(x){(x-filter.bandwidth):(x+filter.bandwidth)})
  sc.jumpNeigh = sort(unique(sc.jumpNeigh))

  assert_that(length(sc.jumpNeigh) < n)

  custom.func <- function(trial){
    set.seed(10*trial)
    y2 = staircase + residual[sample(1:n)]

    #run fused lasso on the new signal
    flasso2 = fusedlasso1d(y2)
    fit2 = coef(flasso2, lambda = lambda)$beta

    filter.val = apply.filter(fit2, bandwidth = filter.bandwidth, return.type = "filter")

    #look at the largest filter value that is more than filter.bandwidth away from a true jump
    max(abs(filter.val[-sc.jumpNeigh]))
  }

  threshold.val = foreach(trial = 1:trials) %dopar% custom.func(trial)

  quantile(threshold.val, prob = quant)
}

.arrange.staircase <- function(fit, tol){
  uniq.val = sort(unique(fit))

  idx.list = lapply(uniq.val, function(x){
    which(fit == x)
  })
  names(idx.list) = uniq.val

  #merge uniq values that less than tol apart
  for(i in (length(uniq.val)-1):1){
    if(abs(uniq.val[i] - uniq.val[i+1]) <= tol) {
      uniq.val = uniq.val[-(i+1)]
      idx.list[[i]] = c(idx.list[[i]], idx.list[[i+1]])
      idx.list[[i+1]] = NULL
    }
  }

  assert_that(all(as.numeric(names(idx.list)) == uniq.val))

  len.vec = sapply(idx.list, length)
  rep(uniq.val, times = len.vec)
}


