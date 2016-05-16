
.plot.helper <- function(jump.location, jump.mean, n, col = "black", 
 lwd = 3, lty = 1, plot.vertical = T){
  assert_that(length(jump.location) == length(jump.mean))
  assert_that(max(jump.location) < n)
  assert_that(all(jump.location == sort(jump.location, decreasing = F)))
  assert_that(jump.location[1] == 1)  

  len = length(jump.location)
  if(len>1){
    for(i in 1:(len-1)){
      lines(x = c(jump.location[i], jump.location[i+1]), y = rep(jump.mean[i],2), 
       col = col, lwd = lwd, lty = lty)
    }
  } 

  #plot the vertical lines
  if(plot.vertical){
    for(i in 2:len){
      lines(x = rep(jump.location[i], 2), y = c(jump.mean[i-1], jump.mean[i]),
       col = col, lwd = lwd, lty = lty)
    }
  }

  #plot the last segment
  lines(x = c(jump.location[len],n), y = rep(jump.mean[len], 2), 
   col = col, lwd = lwd, lty = lty)
  
  invisible()
}


plotfused <- function(jump.mean, jump.location, y, fit, truth = NA, lambda, 
                      num.est.jumps = NA, mse = NA, tol = 1e-7,
                      filter.bandwidth = NA, plotDual = T){
  
  n = length(fit)
  true.seq = form.truth(jump.mean, jump.location, n)
  if(is.na(mse)) mse = compute.mse(fit, true.seq = true.seq)
  if(is.na(num.est.jumps)) num.est.jumps = length(enumerate.jumps(fit))
 
  par(mfrow=c(2,1),mar=c(1,1,1,1))
  
  plot(y, col=rgb(.5,.5,.5), pch=16, cex=1.25)

  .plot.primal(jump.mean, jump.location, y, fit, tol)
  
  #plot dual
  if(plotDual) {
    .plot.dual(jump.location, y, fit, lambda, num.est.jumps, mse, tol)
  } else {
    if(is.na(filter.bandwidth)) filter.bandwidth = ceiling(0.25*log(n)^2)
    
    #WARNING: make truebeta an option
    .plot.filter(fit, filter.bandwidth, jump.mean, jump.location, lambda, mse, 
                 num.est.jumps)
  }
  
  invisible()
}

.extract.location <- function(jump.location, sequence){
  jump.location2 = sapply(jump.location, function(x){max(min(which(
    sequence>=x)),1)-1})
  jump.location2[1] = 1
  jump.location2 = sort(jump.location2)
  
  jump.location2
}

.plot.primal <- function(jump.mean = NA, jump.location = NA, y, 
 fit, tol, truth = NA, verbose = T){
  n = length(y)
  
  if(is.na(truth)){
    true.seq = form.truth(jump.mean, jump.location, n)
  
    #plot truth
    tmp = seq(0, 1,length.out = n)
    jump.location2 = .extract.location(jump.location, tmp)
    .plot.helper(jump.location2, jump.mean, n, col = 2)
  } else {
    truth.split =  split.signal(truth)
    .plot.helper(truth.split$location, truth.split$mean, n, col = 2)
  }

  res.split = split.signal(fit)
  .plot.helper(res.split$location, res.split$mean, n, col = 3)
  
  if(verbose){
    #put text up for a pseudo-y-axis
    tmp.up = round(median(y)+diff(range(y))*.4,2)
    tmp.down = round(median(y)-diff(range(y))*.4,2)
    text(x=0, y=tmp.up, labels=as.character(tmp.up),col = 2)
    text(x=n, y=0, labels=as.character(0),col = 2)
    text(x=0, y=tmp.down, labels=as.character(tmp.down),col = 2)
  }

  invisible()
}

.plot.dual <- function(jump.location, y, fit, lambda, 
                       num.est.jumps, mse, tol){
  tmp = fit-y
  z = cumsum(tmp)
  n = length(y)

  plot(z,ylim=c(-1.5*lambda,1.5*lambda),col=rgb(.5,.5,.5),pch=16)
  lines(z,col="blue",lwd=2)
  lines(x=c(-100,n+100),y=rep(lambda,2),lty=2,lwd=2,col="red")
  lines(x=c(-100,n+100),y=rep(-lambda,2),lty=2,lwd=2,col="red")
  
  tmp = seq(0, 1,length.out = n)
  jump.location2 = .extract.location(jump.location, tmp)
  
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
  text(x=0,y=-1.2*lambda,labels=paste("MSE: ", round(mse,3)," // Lambda: ",round(lambda,2),
                                      " // num.est.jumps: ",num.est.jumps, sep=""),pos=4)
  
  invisible()
  
}

#WARNING: Fix the default for num.est.jumps
.plot.filter <- function(fit, filter.bandwidth, jump.mean, jump.location, 
                         lambda, mse, num.est.jumps = NA, truebeta = NA,
                         verbose = F){
  n = length(fit)
  
  min.dif = min(abs(diff(jump.mean)))
  
  z = apply.filter(fit, filter.bandwidth, min.dif/2, return.type = "filter")
  
  if(all(!is.na(truebeta))){
    truez = apply.filter(truebeta, filter.bandwidth, min.dif/2, 
                         return.type = "filter")
    ylim = c(0, max(min.dif, z, truez))
    plot(truez, ylim = ylim, col = "green", pch = 16)
    
    points(z, col = "blue", pch = 16)
    
  } else {
    ylim = c(0, max(min.dif,z))
    plot(z, ylim = ylim, col = "blue", pch = 16)
  }
  
  
  lines(x = c(-n, 2*n), y = rep(-min.dif/2, 2), lty = 2, lwd = 2, col = "red")
  lines(x = c(-n, 2*n), y = rep(min.dif/2, 2), lty = 2, lwd = 2, col = "red")
  lines(x = c(-n, 2*n), y = rep(0, 2), lty = 2, lwd = 2, col = "red")
  
  n = length(fit)
  #plot the filtered jump locations
  jump.filter = apply.filter(fit, filter.bandwidth, min.dif/2, return.type = "location")
  max.bound = max(abs(z))
  
  for(i in 1:length(jump.filter)){
    lines(x = rep(jump.filter[i],2), y = c(-5*max.bound, 5*max.bound))
  }
  
  
  #plot the true jump locations
  tmp = seq(0, 1,length.out = n)
  jump.location2 = .extract.location(jump.location, tmp)
  
  for(i in 1:length(jump.location2)){
    lines(x = rep(jump.location2[i],2),y = c(-5*max.bound, 5*max.bound), lty = 2, 
          lwd = 2, col = "red")
  }
  
  #put text up for a pseudo-y-axis
  tmp.up = round(median(z)+diff(range(z))*.4, 2)
  tmp.down = round(median(z)-diff(range(z))*.4, 2)
  text(x=0, y=tmp.up, labels=as.character(tmp.up),col="red")
  text(x=0, y=tmp.down, labels=as.character(tmp.down),col="red")
  
  #some basic text on the bottom (mse, lambda, n)
  if(verbose){
    text(x = n, y = 0.7*min(z), labels = paste("MSE: ", round(mse,3), "\nLambda: ", round(lambda,2),
     "\nnum.est.jumps: ", num.est.jumps, 
     "\nFilter Width: ", filter.bandwidth, sep = ""), pos = 2, cex = 0.8)
  }

  invisible()
  
}
