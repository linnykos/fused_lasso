apply.filter <- function(fit, bandwidth, threshold = 0, y = NA, return.type = c("location", "filter")){
  n = length(fit)
  return.type = return.type[1]
  
  jump.loc = enumerate.jumps(fit)

  z =  .compute.filter(fit, bandwidth)

  if(return.type == "filter") return(z)
  
  jump.remain = .compute.candidatejumps(threshold, z, jump.loc, bandwidth, n)
  if(return.type == "location") return(jump.remain)
}

.compute.filter <- function(fit, bandwidth){
  n = length(fit)

  z = sapply(bandwidth:(n-bandwidth), function(x){
    abs(mean(fit[(x+1):(x+bandwidth)]) - mean(fit[(x-bandwidth+1):(x)]))
  })

  #pad the z's
  z = c(rep(0, bandwidth-1), z, rep(0, bandwidth))

  z
}

.compute.candidatejumps <- function(threshold, filter.values, jump.location, bandwidth, n){
  jump.loc2 = c(jump.location, jump.location-bandwidth, jump.location+bandwidth)
  jump.loc2 = jump.loc2[which(jump.loc2 >= bandwidth)]
  jump.loc2 = jump.loc2[which(jump.loc2 <= n - bandwidth)]
  jump.loc2 = c(bandwidth, jump.loc2, n-bandwidth)
  jump.loc2 = sort(unique(jump.loc2))

  filtered.idx = which(abs(filter.values) >= threshold)

  intersect(jump.loc2, filtered.idx)
}

oracle.threshold <- function(fit, bandwidth, truth, threshold.seq = seq(0, 2, length.out = 101)){
  n = length(fit)
  z = .compute.filter(fit, bandwidth)

  jump.loc = enumerate.jumps(fit)
  true.jumps = enumerate.jumps(truth)

  haus.vec = sapply(threshold.seq, function(x){
    jumps.oracle = .compute.candidatejumps(x, z, jump.loc, bandwidth, n)
    compute.hausdorff(true.jumps, jumps.oracle)
  })

  #select the threshold corresponding the smallest haus dist
  #if there are ties, which.min selects the smallest index (i.e., smallest thres)
  threshold.seq[which.min(haus.vec)]
}

classification.quality <- function(true.jumps, estimate.jumps, bandwidth, n){
  #create the augmented set of the true.jumps
  aug.set = function(set){
   unique(as.numeric(unlist(sapply(set, function(x){
    c(max(1, x-bandwidth):min(n, x+bandwidth))
   }))))
  }

  true.neg = setdiff(1:n, aug.set(true.jumps)) #the set of points not close to a true changepoints
  est.pos = aug.set(estimate.jumps) #the set of points close to an estimate changepoint

  #see how many estimated jumps are in the aug.set
  false.pos.ind = any(true.neg %in% estimate.jumps) #is there a true negative in the set of est. jumps
  true.pos.ind = all(true.jumps %in% est.pos) #are all true jumps accounted for

  list(true.pos = true.pos.ind, false.pos = false.pos.ind)
}
