.lower.pos.int = function(x){
  m = length(x)
  splus = sign(x[1])
  zplus = x * splus
  for (i in 1:m){
    zplus[i] = max(min(zplus[i-1],zplus[i]),0)
  }
  
  list(splus,zplus)
}

.lower.int = function(x){
  for.int = .lower.pos.int(x)
  rev.int = .lower.pos.int(rev(x))
  rev.int[[2]] = rev(rev.int[[2]])
  int = for.int[[1]] * for.int[[2]]
  rev.dom = (rev.int[[2]] > for.int[[2]])
  int[rev.dom] = rev.int[[1]]*rev.int[[2]][rev.dom]

  int
}

lower.interpolant <- function(fit, truth){
  jumps = enumerate.jumps(truth, include.endpoints = T)
  z = rep(0, length(fit))

  #compute the interpolant
  for(i in 1:(length(jumps)-1)){

    #first remove the mean
    segment = fit[jumps[i]:(jumps[i]-1)]
    demeaned = segment - mean(segment)

    z[jumps[i]:(jumps[i]-1)] = .lower.int(segment)
  }

  z
}

