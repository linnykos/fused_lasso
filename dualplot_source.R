library(genlasso)
library(assertthat)

generate.problem <- function(n, jump.location, jump.mean, sigma = 0, random.seq = NA){
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



plot.helper <- function(jump.location, jump.mean, n, col = "black", lwd = 3, lty = 1){
  
  len = length(jump.location)
  for(i in 1:(len-1)){
    lines(x=c(jump.location[i],jump.location[i+1]),y=rep(jump.mean[i],2),col=col,lwd=lwd,lty=lty)
  }
  lines(x=c(jump.location[len],n),y=rep(jump.mean[len],2),col=col,lwd=lwd,lty=lty)
  
  invisible()
}


plotfused <- function(jump.mean, jump.location, y, res, lambda, 
                      num.est.jumps, mse = NA, tol = 1e-4){
  par(mfrow=c(2,1),mar=c(1,1,1,1))
  plot(y,col=rgb(.5,.5,.5),pch=16,cex=1.25)
  n = length(y)
  if(is.na(mse)) mse = compute.mse(jump.mean, jump.location, res)
  
  #plot truth
  tmp = seq(0,1,length.out=n)
  jump.location2 = sapply(jump.location,function(x){max(min(which(tmp>=x)),1)-1})
  jump.location2[1] = 1
  jump.location2 = sort(jump.location2)
  plot.helper(jump.location2, jump.mean, n, col="red")
  
  
  tmpdiff = diff(res)
  tmpdiff[abs(tmpdiff)<tol] = 0
  res.jumploc = c(1,(which(tmpdiff != 0)))
  res.jumpmean = res[res.jumploc+1]
  res.jumploc = sort(res.jumploc)
  plot.helper(res.jumploc, res.jumpmean, n, col=rgb(0,1,0))
  
  #put text up for a pseudo-y-axis
  tmp.up = round(median(y)+diff(range(y))*.4,2)
  tmp.down = round(median(y)-diff(range(y))*.4,2)
  text(x=0, y=tmp.up, labels=as.character(tmp.up),col="red")
  text(x=n, y=0, labels=as.character(0),col="red")
  text(x=0, y=tmp.down, labels=as.character(tmp.down),col="red")
  
  #plot dual
  tmp = res-y
  z = cumsum(tmp)
  plot(z,ylim=c(-1.5*lambda,1.5*lambda),col=rgb(.5,.5,.5),pch=16)
  lines(z,col="blue",lwd=2)
  lines(x=c(-100,n+100),y=rep(lambda,2),lty=2,lwd=2,col="red")
  lines(x=c(-100,n+100),y=rep(-lambda,2),lty=2,lwd=2,col="red")
  for(i in 1:length(jump.location2)){
    lines(x=rep(jump.location2[i],2),y=c(-5*lambda,5*lambda),lty=2,lwd=2,col="red")
  }
  
  tmp = which(abs(abs(z)-lambda)<tol)
  for(i in 1:length(tmp)){
    if(z[tmp[i]]>0) {
      lines(x=rep(tmp[i],2),y=c(lambda,5*lambda))
    } else {
      lines(x=rep(tmp[i],2),y=c(-lambda,-5*lambda))
    } 
  }
  
  #some basic text on the bottom (mse, lambda, n)
  text(x=0,y=-1.2*lambda,labels=paste("MSE: ", round(mse,3)," // Lambda: ",round(lambda,2)," // num.est.jumps: ",num.est.jumps,sep=""),pos=4)
  
  invisible()
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
    num.est.jumps = count.jumps(res)
  } else {
    tmp = coef(fres,df=num.est.jumps)
    res = tmp$beta
    lambda = tmp$lambda
  }
  mse = compute.mse(jump.mean, jump.location, res)
  
  if(plot.res) plotfused(jump.mean,jump.location,y,res,lambda,num.est.jumps,mse)

  ret = list()
  ret[["coef"]] = res
  ret[["lambda"]] = lambda
  ret[["jumps"]] = num.est.jumps
  ret[["mse"]] = mse
  if(return.fusedlasso) ret[["fusedlasso"]] = fres
  
  return(ret)
}

compute.mse <- function(jump.mean, jump.location, res){
  n = length(res)
  
  tmp = seq(0,1,length.out=n)
  jump.location2 = sapply(jump.location,function(x){max(min(which(tmp>=x)),1)-1})
  jump.location2[1] = 1
  jump.location2 = sort(jump.location2)
  jump.location3 = c(jump.location2,n)
  jump.location3[1] = 0

  true.seq = rep(jump.mean,times=diff(jump.location3))
  
  sum((true.seq-res)^2)/n
}

count.jumps <- function(res, tol=1e-3){
  length(which(abs(diff(res))>tol))
}

#dummy function temporarily
kkt_checker <- function(){
  
  tmp = res-y
  n = length(y)
  z = cumsum(tmp)
  idx = which(abs(abs(z)-lambda)<tol)
  tmpdiff = diff(res)
  tmpdiff[abs(tmpdiff)<tol] = 0
  res.jumploc = c(1,(which(tmpdiff != 0)+1))
  res.jumpmean = res[res.jumploc]
  svec = sign(diff(res.jumpmean))
  #NEED TO BUMP BY ONE INDEX (z should be a vector of length n+1)
  svec[10]=1
  tmpmean = rep(0,length(svec)+1)
  tmpmean[1] = 1/idx[1]*sum(y[1:idx[1]])+lambda/idx[1]*svec[1]
  for(i in 2:(length(svec))){
    tmpmean[i] = 1/(idx[i]-idx[i-1])*sum(y[(idx[i-1]+1):idx[i]])+lambda/(idx[i]-idx[i-1])*(svec[i]-svec[i-1])
  }
  i = length(svec)+1
  tmpmean[i] = 1/(n-idx[i-1])*sum(y[(idx[i-1]+1):n])+lambda/(n-idx[i-1])*(-svec[i-1])
  #reconstruct the beta vec
  tmp = diff(c(1,res.jumploc[-1],n+1))
  tmpres = rep(tmpmean,tmp)
  #check the z's
  tmp = tmpres-y
  tmpz = cumsum(tmp)
  #check KKT conditions
  print(paste("Max z: ",max(abs(tmpz)),sep=""))
  for(i in 1:length(svec)){
    print(paste("Check jump ",i,": ",tmpz[idx[i]],"=",lambda*sign(tmpres[idx[i]+1]-tmpres[idx[i]]),sep=""))
  }
  print(paste("Check sum: ",sum(y),"=",sum(tmpres),sep=""))
  print(rbind(tmpmean,c(0,svec)))
  plot(z)
  plot(tmpz)
  lines(x=c(-1000,1000),y=rep(lambda,2),col="red",lty=2,lwd=2)
  lines(x=c(-1000,1000),y=rep(-lambda,2),col="red",lty=2,lwd=2)
  
}