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


plotfused <- function(jump.mean,jump.location,y,res,lambda, tol=1e-3){
  par(mfrow=c(2,1),mar=c(1,1,1,1))
  plot(y,col=rgb(.5,.5,.5),pch=16,cex=1.25)
  n = length(y)
  
  #plot truth
  tmp = seq(0,1,length.out=n)
  jump.location = sapply(jump.location,function(x){max(min(which(tmp>=x)),1)-1})
  jump.location[1] = 1
  jump.location = sort(jump.location)
  plot.helper(jump.location, jump.mean, n, col="red")
  
  
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
  text(x=0, y=0, labels=as.character(0),col="red")
  text(x=0, y=tmp.down, labels=as.character(tmp.down),col="red")
  
  #plot dual
  tmp = res-y
  z = cumsum(tmp)
  plot(z,ylim=c(-1.5*lambda,1.5*lambda),col=rgb(.5,.5,.5),pch=16)
  lines(z,col="blue",lwd=2)
  lines(x=c(-100,n+100),y=rep(lambda,2),lty=2,lwd=2,col="red")
  lines(x=c(-100,n+100),y=rep(-lambda,2),lty=2,lwd=2,col="red")
  for(i in 1:length(jump.location)){
    lines(x=rep(jump.location[i],2),y=c(-5*lambda,5*lambda),lty=2,lwd=2,col="red")
  }
  
  tmp = which(abs(abs(z)-lambda)<tol)
  for(i in 1:length(tmp)){
    if(z[tmp[i]]>0) {
      lines(x=rep(tmp[i],2),y=c(lambda,5*lambda))
    } else {
      lines(x=rep(tmp[i],2),y=c(-lambda,-5*lambda))
    } 
  }
  
  invisible()
}

dualplot_suite <- function(y, lambda, jump.mean, jump.location){
  fres = fusedlasso1d(y)
  res = coef(fres,lambda=lambda)$beta
  plotfused(jump.mean,jump.location,y,res,lambda)
  invisible()
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