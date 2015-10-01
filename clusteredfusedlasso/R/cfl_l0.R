cfl.1D.l0 = function(y,  K, lam, niter=50, verbose=FALSE, multiplier = 1, initializer = "l1") {
  
  cl = match.call()
  
  k = K
  lam2 = lam
  n = length(y)
  
  info = matrix(NA,ncol=6,nrow=1)
  colnames(info) = c("n","K","niter","init.lambda","final.lambda","multiplier")
  info[1:6] = c(n,K,niter,lam,lam*multiplier^niter,multiplier)
  attr(info,"initializer") = initializer
  
  obj.vec = rep(NA,niter+1)
  mse.vec = rep(NA,niter+1)
  true.mse.vec = rep(NA,niter+1)
  
  # initialize the groups and beta
  if(initializer=="l1"){
    out = fusedlasso1d(y)
    tmp = coef.genlasso(out,lambda=lam)$beta
  } else {
    out = L0Seg(y, lambda2 = lam)
    tmp = out$fit
  }
  mse.fused = 1/n*sum((y-tmp)^2)
  out2 = Ckmeans.1d.dp(tmp,k=k)
  beta = out2$centers
  g = out2$cluster
  
  obj.vec[1] = crit.1D.l0(y, beta, g, lam)
  mse.vec[1] = mse.1D(y, beta, g, n)
  
  # begin alternating optimization
  for (i in 1:niter) {
    
    if (verbose) cat(paste(i,":  ",sep=""))
    beta = sfl.1D.l0(y,g,n,k,lam)
    #check to see if any clusters disappeared
    k = k - check.unique(beta,k)
    if (verbose) cat(paste(round(crit.1D.l0(y,beta,g,lam2),3)," "))
    g = dfl.1D.l0(y,beta,n,k,lam2)
    if (verbose) cat(paste(round(crit.1D.l0(y,beta,g,lam2),3)," "),fill=TRUE)
    
    obj.vec[i+1] = crit.1D.l0(y, beta, g, lam)
    mse.vec[i+1] = mse.1D(y,beta, g, n)
    
    lam2 = lam2*multiplier
  }
  
  list(call = cl, 
       beta=beta, g=g, 
       obj.vec=obj.vec, mse.vec=mse.vec, 
       info = info,
       y=y, 
       fused.sol = tmp, init.sol = out2$centers[out2$cluster], 
       mse.fused = mse.fused)
}

# this is a wrapper function sets things up in the 1D case
sfl.1D.l0 = function(y, g, n, K, lam) {
  jj = rep(TRUE,K)
  for (j in 1:K) jj[j] = any(g==j)
  
  beta = rep(0,K)
  beta[jj] = sfl.l0(y,cumsum(jj)[g],n, sum(jj),lam)
  
  beta
}

# small fused lasso
sfl.l0 = function(y, g, n, K, lam) {
  if (K==1) {
    return(mean(y))
  }
  
  beta = rep(NA,K)
  for(i in 1:K){
    beta[i] = mean(y[which(g==i)])
  }
  
  beta
}

# discrete fused lasso in 1D
# (the C library "dfl.so" must be loaded before you call this function)
dfl.1D.l0 = function(y, beta, n, K, lam) {
  o = .C("dfl1Dl0", n=as.integer(n), y=as.double(y),
         K=as.integer(K), b=as.double(beta), lam=as.double(lam),
         g=as.integer(numeric(n)),dup=FALSE)
  
  o$g+1
}

# the criterion in 1D
crit.1D.l0 = function(y, beta, g, lam) {
  1/2*sum((y-beta[g])^2) + lam*length(which(abs(diff(beta[g]))>0.001))
}