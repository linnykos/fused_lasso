generate.problem <- function(n,jump.location,jump.mean,sigma, scale=FALSE){
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
#   
#   if(scale){
#     ymax = max(y)
#     ymin = min(y)
#     y = (y-ymin)/(ymax-ymin)
#     C = ymax-ymin
#   } else {
#     C = NA
#   }
#   
#   list(y=y, C=C)
  y
}

plot.helper <- function(jump.location, jump.mean, col = "black", lwd = 3, lty = 1){
  
  len = length(jump.location)
  for(i in 1:(len-1)){
    lines(x=c(jump.location[i],jump.location[i+1]),y=rep(jump.mean[i],2),col=col,lwd=lwd,lty=lty)
  }
  lines(x=c(jump.location[len],n),y=rep(jump.mean[len],2),col=col,lwd=lwd,lty=lty)
  
  invisible()
}

#res1 should be the initial solution from fused lasso
#res2 (optional) should be the new initial solution after running K-means
#res3 (optional) should be the final solution after running alternating minimization
plotspecial <- function(res1, res2=TRUE, res3=TRUE, lambda, 
                        jump.location = NA, jump.mean = NA, tol = 0.001,
                        draw.kmeans = TRUE, onlyfused = FALSE){
  
  vec.coef = res1$beta[res1$g] #coef.genlasso(res1,lambda=1)$beta
  l1 = NA
  l2 = NA
  linf = NA
  if(!is.na(jump.location[1])) {
    evalmat = eval.coef(res, jump.location, jump.mean, sigma)
    l1 = round(evalmat[3,1],4)
    l2 =  round(evalmat[4,1],4)
    l2fused =  round(evalmat[4,2],4)
    l2init =  round(evalmat[4,3],4)
    linf =  round(evalmat[5,1],4)
  }
  
  n = length(res1$g)
  truth = construct.truth(jump.location, jump.mean, n)
  h1 = round(hausdorff.dist(res1$beta[res$g],truth),3)
  h2 =  round(hausdorff.dist(res1$fused.sol,truth),3)
  h3 =  round(hausdorff.dist(res1$init.sol,truth),3)
  
  plot(res1$y,pch=16,col=rgb(.5,.5,.5),xlab="Position",ylab="Y",main=paste("Lambda: ",res$info[4],"// L2-Fused1: ",l2fused," (h: ",h2,")","//\n L2-Clust Fused1: ",l2init," (h: ",h3,")","// L2-cF0Lasso: ",l2," (h: ",h1,")",sep=""))
  
  legend("topleft", c("Points","Truth","f1Lasso","clust. f1Lasso","cfLasso"), bty="n", fill=c("gray","red","green","orange","blue"),cex=0.75)
  
  #plot the truth
  #convert jump.location into integer indices
  tmp = seq(0,1,length.out=n)
  jump.location = sapply(jump.location,function(x){min(which(tmp>=x))})
  jump.location[1] = 1
  plot.helper(jump.location, jump.mean, col="red")
  
  if(!onlyfused){
    #plot the res1
    tmpdiff = diff(vec.coef)
    tmpdiff[abs(tmpdiff)<tol] = 0
    res1.jumploc = c(1,(which(tmpdiff != 0)+1))
    res1.jumpmean = vec.coef[res1.jumploc]
    plot.helper(res1.jumploc, res1.jumpmean, col="blue")
  }
  
  
  if(res2) res2 = as.numeric(res1$fused.sol)
  tmpdiff = diff(res2)
  tmpdiff[abs(tmpdiff)<tol] = 0
  res2.jumploc = c(1,(which(tmpdiff != 0)+1))
  res2.jumpmean = res2[res2.jumploc]
  plot.helper(res2.jumploc, res2.jumpmean, col="green")
  
  if(!onlyfused){
    if(res3) res3 = res1$init.sol
    tmpdiff = diff(res3)
    tmpdiff[abs(tmpdiff)<tol] = 0
    res3.jumploc = c(1,(which(tmpdiff != 0)+1))
    res3.jumpmean = res3[res3.jumploc]
    plot.helper(res3.jumploc, res3.jumpmean, col="orange")
  }

  
  if(draw.kmeans){
    tmp = Ckmeans.1d.dp(res$y,res$info[2])$centers
    for(i in 1:length(tmp)){
      lines(x=c(-100,res$info[1]+100),y=rep(tmp[i],2),col="red",lty="dotted")
    }
  }
  
  invisible() 
}


plotspecial.l0 <- function(res1, res2=NA, lambda, 
                        jump.location = NA, jump.mean = NA, tol = 0.001,
                        draw.kmeans = TRUE, onlyfused = FALSE){
  
  vec.coef = res1$beta[res1$g] #coef.genlasso(res1,lambda=1)$beta
  l1 = NA
  l2 = NA
  linf = NA
  if(!is.na(jump.location[1])) {
    evalmat = eval.coef(res, jump.location, jump.mean, sigma)
    l1 = round(evalmat[3,1],4)
    l2 =  round(evalmat[4,1],4)
    
    evalmat = eval.coef(res2, jump.location, jump.mean, sigma)
    l2fused =  round(evalmat[4,2],4)
    l2init =  round(evalmat[4,3],4)
  }

  n = length(res1$g)
  truth = construct.truth(jump.location, jump.mean, n)
  h1 = round(hausdorff.dist(res1$beta[res$g],truth),3)
  h2 =  round(hausdorff.dist(res2$fused.sol,truth),3)
  h3 =  round(hausdorff.dist(res2$init.sol,truth),3)

  plot(res1$y,pch=16,col=rgb(.5,.5,.5),xlab="Position",ylab="Y",main=paste("Lambda: ",res$info[4],"// L2-Fused0: ",l2fused," (h: ",h2,")","//\n L2-Clust Fused0: ",l2init," (h: ",h3,")","// L2-cF0Lasso: ",l2," (h: ",h1,")",sep=""))
  
  legend("topleft", c("Points","Truth","f0Lasso","clust. f0Lasso","cf0Lasso"), bty="n", fill=c("gray","red","green","orange","blue"),cex=0.75)
  
  #plot the truth
  #convert jump.location into integer indices
  tmp = seq(0,1,length.out=n)
  jump.location = sapply(jump.location,function(x){min(which(tmp>=x))})
  jump.location[1] = 1
  plot.helper(jump.location, jump.mean, col="red")
  
  #plot the res1
  tmpdiff = diff(vec.coef)
  tmpdiff[abs(tmpdiff)<tol] = 0
  res1.jumploc = c(1,(which(tmpdiff != 0)+1))
  res1.jumpmean = vec.coef[res1.jumploc]
  plot.helper(res1.jumploc, res1.jumpmean, col="blue")  
  
  res2b = res2$fused.sol
  tmpdiff = diff(res2b)
  tmpdiff[abs(tmpdiff)<tol] = 0
  res2.jumploc = c(1,(which(tmpdiff != 0)+1))
  res2.jumpmean = res2b[res2.jumploc]
  plot.helper(res2.jumploc, res2.jumpmean, col="green")
  
  res3b = res2$init.sol
  tmpdiff = diff(res3b)
  tmpdiff[abs(tmpdiff)<tol] = 0
  res3.jumploc = c(1,(which(tmpdiff != 0)+1))
  res3.jumpmean = res3b[res3.jumploc]
  plot.helper(res3.jumploc, res3.jumpmean, col="orange")

  
  
  if(draw.kmeans){
    tmp = Ckmeans.1d.dp(res$y,res$info[2])$centers
    for(i in 1:length(tmp)){
      lines(x=c(-100,res$info[1]+100),y=rep(tmp[i],2),col="red",lty="dotted")
    }
  }
  
  invisible() 
}

#the "test set" which draws new data from same distribution and evaluates
eval.coef <- function(res, jump.location, jump.mean, sigma, lam = NA){
  n = length(res$y)

  #reconstruct truth
  tmp = seq(0,1,length.out=n)
  jump.location = sapply(jump.location,function(x){min(which(tmp>=x))})
  jump.location[1] = 1
  truth = rep(0,n)
  for(i in 1:(length(jump.location)-1)){
    tmp = jump.location[i]:(jump.location[i+1]-1)
    truth[tmp] = jump.mean[i]
  }
  tmp = jump.location[length(jump.location)]:n
  truth[tmp] = jump.mean[length(jump.mean)]
  
  mat = matrix(NA,ncol=4,nrow=5)
  colnames(mat) = c("Final","Fused", "Initial","Truth")
  rownames(mat) = c("MSE","Obj","L1-Truth","L2-Truth","Linf-Truth")

  mat[1,1] = res$mse.vec[length(res$mse.vec)]
  mat[1,2] = res$mse.fused
  mat[1,3] = res$mse.vec[1]
  mat[1,4] = 1/n*sum((res$y-truth)^2)
  
  
  if(is.na(lam)) lam = res$info[5]
  #compute objective functions
  y = res$y
  mat[2,1] = 1/2*sum((y-res$beta[res$g])^2) + lam*length(which(abs(diff(res$beta[res$g]))>0.001))
  mat[2,2] = 1/2*sum((y-res$fused.sol)^2) + lam*length(which(abs(diff(res$fused.sol))>0.001))
  mat[2,3] = 1/2*sum((y-res$init.sol)^2) + lam*length(which(abs(diff(res$init.sol))>0.001))
  mat[2,4] = 1/2*sum((y-truth)^2) + lam*length(which(abs(diff(truth))>0.001))
  
  #compute the l1 distance from true mean
  mat[3,1] = 1/n*sum(abs(res$beta[res$g]-truth))
  mat[3,2] = 1/n*sum(abs(res$fused.sol-truth))
  mat[3,3] = 1/n*sum(abs(res$init.sol-truth))
  
  #compute the l2 distance from true mean
  mat[4,1] = 1/n*sum((res$beta[res$g]-truth)^2)
  mat[4,2] = 1/n*sum((res$fused.sol-truth)^2)
  mat[4,3] = 1/n*sum((res$init.sol-truth)^2)
  
  #compute the l_inf distance from true mean
  mat[5,1] = max(abs(res$beta[res$g]-truth))
  mat[5,2] = max(abs(res$fused.sol-truth))
  mat[5,3] = max(abs(res$init.sol-truth))

  mat
  
}


construct.truth <- function(jump.location, jump.mean, n){
  tmp = seq(0,1,length.out=n)
  jump.location = sapply(jump.location,function(x){min(which(tmp>=x))})
  jump.location[1] = 1
  truth = rep(0,n)
  for(i in 1:(length(jump.location)-1)){
    tmp = jump.location[i]:(jump.location[i+1]-1)
    truth[tmp] = jump.mean[i]
  }
  tmp = jump.location[length(jump.location)]:n
  truth[tmp] = jump.mean[length(jump.mean)]
  
  truth
}

hausdorff.dist <- function(res1, res2, tol = 0.001){
  n = length(res1)
  
  tmpdiff = diff(res1)
  tmpdiff[abs(tmpdiff)<tol] = 0
  res1.jumploc = c(1,(which(tmpdiff != 0)+1))
  
  tmpdiff = diff(res2)
  tmpdiff[abs(tmpdiff)<tol] = 0
  res2.jumploc = c(1,(which(tmpdiff != 0)+1))

  n1 = length(res1.jumploc)
  n2 = length(res2.jumploc)
  
  vec1 = sapply(res1.jumploc,function(x){min(abs(res2.jumploc-x))})
  vec2 = sapply(res2.jumploc,function(x){min(abs(res1.jumploc-x))})
  
  max(c(vec1,vec2))
}