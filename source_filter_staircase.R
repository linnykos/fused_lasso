.staircasecontrol <- setClass("staircasecontrol", representation(trials = 
 "numeric", quant = "numeric", tol = "numeric", type = "character"),
 prototype(trials = 100, quant = 0.95, tol = 1e-6, type = "flip"))

#rearrange the fit to be a big staircase
staircase.threshold <- function(y, fit, filter.bandwidth, lambda, controls = 
 list(trials = 50, quant = 0.95, tol = 1e-6, type = "flip")){

  con = .convert.list2control(controls, "staircasecontrol")

  assert_that(length(y) == length(fit))

  n = length(y)
  residual = y - fit

  if(con@type == "rearrange") {
    baseline = .arrange.staircase(fit)
  } else if(con@type == "flip") {
    baseline = .arrange.flip(fit, con@tol)
  } else if(con@type == "original"){
    baseline = fit
  }

  sc.jumps = enumerate.jumps(baseline, con@tol)

  sc.jumpNeigh = sapply(sc.jumps, function(x){
   max(1,(x-filter.bandwidth)) : min((x+filter.bandwidth), n)
  })

  sc.jumpNeigh = sort(unique(as.numeric(unlist(sc.jumpNeigh))))

  if(length(sc.jumpNeigh) >= n) return(rep(NA, length(con@quant)))

  custom.func <- function(trial){
    set.seed(10*trial)
    y2 = baseline + residual[sample(1:n)]

    #run fused lasso on the new signal
    flasso2 = fusedlasso1d(y2)
    fit2 = coef(flasso2, lambda = lambda)$beta

    filter.val = apply.filter(fit2, bandwidth = filter.bandwidth, 
     return.type = "filter")

    #look at the largest filter value that is more than filter.bandwidth away from a true jump
    max(abs(filter.val[-sc.jumpNeigh]))
  }

  threshold.val = foreach(trial = 1:con@trials) %dopar% custom.func(trial)

  quantile(unlist(threshold.val), prob = con@quant)
}

.arrange.staircase <- function(fit){
  uniq.val = sort(unique(fit))

  idx.list = lapply(uniq.val, function(x){
    which(fit == x)
  })
 
  len.vec = sapply(idx.list, length)
  rep(uniq.val, times = len.vec)
}

.arrange.flip <- function(fit, tol = 0){
  jumps = enumerate.jumps(fit, tol)
  jumps = c(1, jumps, length(fit))
  newfit = fit

  #on iteration i, figure out the difference between i-1 and i
  #and replace all the values at i with abs val
  for(i in 2:(length(jumps)-1)){
    diff.val = abs(fit[jumps[i]] - fit[jumps[i-1]])
    len = jumps[i+1] - jumps[i]

    newfit[jumps[i]:(jumps[i+1]-1)] = rep(newfit[jumps[i-1]]+diff.val, len)
  }

  newfit
}
