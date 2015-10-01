# clustered fused lasso in 1D
cfl.1D = function(y,  K, lam, niter=50, verbose=FALSE, multiplier = 1) {
  
  cl = match.call()
  
  k = K
  lam2 = lam
  n = length(y)
  
  info = matrix(NA,ncol=6,nrow=1)
  colnames(info) = c("n","K","niter","init.lambda","final.lambda","multiplier")
  info[1:6] = c(n,K,niter,lam,lam*multiplier^niter,multiplier)
  
  obj.vec = rep(NA,niter+1)
  mse.vec = rep(NA,niter+1)
  true.mse.vec = rep(NA,niter+1)

  # initialize the groups and beta
  out = fusedlasso1d(y)
  tmp = coef.genlasso(out,lambda=lam)$beta
  mse.fused = 1/n*sum((y-tmp)^2)
  out2 = Ckmeans.1d.dp(tmp,k=k)
  beta = out2$centers
  g = out2$cluster
  
  obj.vec[1] = crit.1D(y, beta, g, lam)
  mse.vec[1] = mse.1D(y, beta, g, n)

  # begin alternating optimization
  for (i in 1:niter) {
    
    if (verbose) cat(paste(i,":  ",sep=""))
    beta = sfl.1D(y,g,n,k,lam)
    #check to see if any clusters disappeared
    k = k - check.unique(beta,k)
    if (verbose) cat(paste(round(crit.1D(y,beta,g,lam2),3)," "))
    g = dfl.1D(y,beta,n,k,lam2)
    if (verbose) cat(paste(round(crit.1D(y,beta,g,lam2),3)," "),fill=TRUE)
    
    obj.vec[i+1] = crit.1D(y, beta, g, lam)
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
sfl.1D = function(y, g, n, K, lam) {
  jj = rep(TRUE,K)
  for (j in 1:K) jj[j] = any(g==j)
  
  beta = rep(0,K)
  beta[jj] = sfl(y,cumsum(jj)[g],n, sum(jj),lam)
  
  beta
}

# small fused lasso
# this is a QP with K*(K-1)/2 variables and K*(K-1) box constraints
sfl = function(y, g, n, K, lam) {
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
   cc[cc==0] = 0.000001
  
  D = Dtmp%*%diag(1/sqrt(nn))

  #alg fix
  Qtmp = D%*%t(D)
  Qtmp = Qtmp+diag(K2)*0.01
  dtmp = -D%*%z
  Atmp = t(rbind(diag(K2),-diag(K2)))
  btmp = c(-lam*cc,-lam*cc)
  
  uvec = solve.QP(Dmat = Qtmp, dvec = dtmp, Amat = Atmp, bvec = btmp)$solution
  
  alpha = z + t(D)%*%uvec
  
  beta = alpha/sqrt(nn)
  
  beta
}

# discrete fused lasso in 1D
# (the C library "dfl.so" must be loaded before you call this function)
dfl.1D = function(y, beta, n, K, lam) {
  o = .C("dfl1D", n=as.integer(n), y=as.double(y),
    K=as.integer(K), b=as.double(beta), lam=as.double(lam),
    g=as.integer(numeric(n)),dup=FALSE)
  
  o$g+1
}

# the criterion in 1D
crit.1D = function(y, beta, g, lam) {
  1/2*sum((y-beta[g])^2) + lam*sum(abs(diff(beta[g])))
}

# the criterion in 1D
mse.1D = function(y, beta, g, n) {
  1/n*sum((y-beta[g])^2)
}

check.unique = function(beta,k,tol=0.001){
  tmp = diff(sort(beta))
  tmp2 = which(abs(tmp)<tol)
  if(length(tmp2)>0) print(paste("CLUSTER WENT DOWN BY ",length(tmp2),sep=""))
  
  length(tmp2)
}