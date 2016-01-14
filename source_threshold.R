threshold.fit <- function(fit, num.jumps){
  num.est.jumps = count.jumps(fit)
  if(num.jumps >= num.est.jumps) return(fit)
  
  dif = abs(diff(fit))
  idx = c(1,sort(order(dif,decreasing=TRUE)[1:num.jumps]),length(y)+1)
  len = length(idx)
  
  #refit the function
  #shrink the beta vector to be smaller
  fit.avg = unlist(sapply(1:(len-1),function(x){
    rep(mean(fit[idx[x]:idx[x+1]-1]), idx[x+1]-idx[x])
  }))
  
  fit.avg
}

#based on the num.jumps largest jumps in f, return the average y value
# between the largest jumps and the indices
threshold.fit.to.y <- function(fit, y, num.jumps){
  num.est.jumps = count.jumps(fit)
  if(num.jumps >= num.est.jumps) return(fit)
  
  dif = abs(diff(fit))
  idx = c(1,sort(order(dif,decreasing=TRUE)[1:num.jumps]),length(y)+1)
  len = length(idx)
  
  #refit the function
  #shrink the beta vector to be smaller
  y.avg = unlist(sapply(1:(len-1),function(x){
    mean(y[idx[x]:idx[x+1]-1])
  }))
  
  list(y.avg = y.avg, idx = idx)
}

refit.estimate <- function(fit, y, num.jumps, lambda=1){
  assert_that(num.jumps>0)
  
  thres = threshold.fit.to.y(fit, y, num.jumps)
  idx = thres$idx
  idx[1] = 0
  idx[length(idx)] = length(y)
  
  idx2 = rep(1:(num.jumps+1), times=diff(idx))
  beta = sfl.1D(y, idx2, lambda)
  
  #reformat the coefficients
  beta.long = rep(beta, times=diff(idx))
  beta.long
}


# this is a wrapper function sets things up in the 1D case
sfl.1D = function(y, g, lam, n = length(y), K=length(unique(g))) {
  jj = rep(TRUE,K)
  for (j in 1:K) jj[j] = any(g==j)
  
  beta = rep(0,K)
  beta[jj] = sfl(y,cumsum(jj)[g],n, sum(jj),lam)
  
  beta
}

# small fused lasso
# this is a QP with K*(K-1)/2 variables and K*(K-1) box constraints
sfl = function(y, g, n, K, lam, tol=1e-5) {
  if (K==1) {
    return(mean(y))
  }
  
  nn = as.numeric(table(g))
  yy = rep(0,K)
  for (j in 1:K) {
    yy[j] = mean(y[g==j])
  }
  
  z = sqrt(nn)*yy
  
  K2 = K*(K-1)/2
  Dtmp = matrix(0,nrow=K2,ncol=K)
  i = 1
  for (j in 1:(K-1)) {
    for (l in (j+1):K) {
      Dtmp[i,j] = -1
      Dtmp[i,l] = 1
      i = i+1
    }
  }
  
  cc = rep(0,K2)
  for (i in 1:(n-1)) {
    if (g[i]!=g[i+1]) {
      j = min(g[i],g[i+1])
      l = max(g[i],g[i+1])
      ii = (j-1)*K-(j-1)*j/2 + (l-j)
      cc[ii] = cc[ii]+1
    }
  }
  
  #   ##ALG FIX: 
  cc[cc==0] = tol
  
  D = Dtmp%*%diag(1/sqrt(nn))
  
  #alg fix
  Qtmp = D%*%t(D)
  Qtmp = Qtmp+diag(K2)*0.01
  dtmp = -D%*%z
  Atmp = t(rbind(diag(K2),-diag(K2)))
  btmp = c(-lam*cc,-lam*cc)
  Qtmp = as.matrix(Qtmp, "sparseMatrix")
  Atmp = as.matrix(Atmp, "sparseMatrix")
  
  uvec = solve.QP(Dmat = Qtmp, dvec = dtmp, Amat = Atmp, bvec = btmp)$solution
  
  alpha = z + t(D)%*%uvec
  
  beta = alpha/sqrt(nn)
  
  beta
}
