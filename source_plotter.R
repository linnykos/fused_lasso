
plot.helper <- function(jump.location, jump.mean, n, col = "black", lwd = 3, lty = 1){
  
  len = length(jump.location)
  if(len>1){
    for(i in 1:(len-1)){
      lines(x=c(jump.location[i],jump.location[i+1]),y=rep(jump.mean[i],2),col=col,lwd=lwd,lty=lty)
    }
  } 
  lines(x=c(jump.location[len],n),y=rep(jump.mean[len],2),col=col,lwd=lwd,lty=lty)
  
  invisible()
}


plotfused <- function(jump.mean, jump.location, y, res, lambda, 
                      num.est.jumps, mse = NA, tol = 1e-7){
  par(mfrow=c(2,1),mar=c(1,1,1,1))
  plot(y,col=rgb(.5,.5,.5),pch=16,cex=1.25)
  n = length(y)
  
  true.seq = form.truth(jump.mean, jump.location, n)
  if(is.na(mse)) mse = compute.mse(res, true.seq = true.seq)
  
  #plot truth
  tmp = seq(0,1,length.out=n)
  jump.location2 = sapply(jump.location,function(x){max(min(which(tmp>=x)),1)-1})
  jump.location2[1] = 1
  jump.location2 = sort(jump.location2)
  plot.helper(jump.location2, jump.mean, n, col="red")
  
  
  tmpdiff = diff(res)
  tmpdiff[abs(tmpdiff)<tol] = 0
  res.jumploc = c(1,(which(abs(tmpdiff)>tol)))
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
