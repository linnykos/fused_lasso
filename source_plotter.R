
.plot.helper <- function(jump.location, jump.mean, n, col = "black", 
 lwd = 3, lty = 1, plot.vertical = T){
  assert_that(length(jump.location) == length(jump.mean))
  assert_that(max(jump.location) < n)
  assert_that(all(jump.location == sort(jump.location, decreasing = F)))
  assert_that(jump.location[1] == 1)  

  len = length(jump.location)
  if(len>1){
    for(i in 1:(len-1)){
      lines(x = c(jump.location[i], jump.location[i+1]), 
       y = rep(jump.mean[i],2), col = col, lwd = lwd, lty = lty)
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
                      filter.bandwidth = NA, threshold = NA, plotDual = T){
  
  n = length(fit)
  true.seq = form.truth(jump.mean, jump.location, n)
  if(is.na(mse)) mse = compute.mse(fit, true.seq = true.seq)
  if(is.na(num.est.jumps)) num.est.jumps = length(enumerate.jumps(fit))
 
  par(mfrow=c(2,1),mar=c(2,4,1,1))
  
  plot(y, col=rgb(.7,.7,.7), pch=16, cex=1.25, ylab = "Data values", 
   cex.axis = .8, cex.lab = .8)

  .plot.primal(jump.mean, jump.location, y, fit, tol)
  
  #plot dual
  if(plotDual) {
    .plot.dual(jump.location, y, fit, lambda, num.est.jumps, mse, tol)
  } else {
    if(is.na(filter.bandwidth)) filter.bandwidth = ceiling(0.25*log(n)^2)
    
    #WARNING: make truebeta an option
    .plot.filter(fit, filter.bandwidth, jump.mean, jump.location, lambda, mse, 
      threshold, num.est.jumps)
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
 fit, tol, truth = NA, verbose = F){
  n = length(y)
  
  if(all(is.na(truth))){
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
  .plot.helper(res.split$location, res.split$mean, n, col = "blue2")
  
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
 num.est.jumps, mse, tol, verbose = F){
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
  if(verobse){
    text(x=0,y=-1.2*lambda,labels=paste("MSE: ", round(mse,3),
     " // Lambda: ",round(lambda,2), " // num.est.jumps: ",
     num.est.jumps, sep=""),pos=4)
  }  

  invisible()
  
}

#WARNING: Fix the default for num.est.jumps
.plot.filter <- function(fit, filter.bandwidth, jump.mean, jump.location, 
 lambda, mse, threshold = NA, num.est.jumps = NA, truebeta = NA, verbose = F){

  n = length(fit)
  
  if(is.na(threshold)) threshold = min(abs(diff(jump.mean)))/2
  
  z = apply.filter(fit, filter.bandwidth, threshold, return.type = "filter")
 
  #plot the filter values of the true beta
  if(all(!is.na(truebeta))){
    truez = apply.filter(truebeta, filter.bandwidth, threshold, 
     return.type = "filter")
    ylim = c(0, max(threshold, z, truez))
    plot(truez, ylim = ylim, col = "blue2", pch = 16, ylab = "Filter values", 
     cex.axis = .8, cex.lab = .8)
    
    points(z, col = "blue2", pch = 16)
    
  } else {
    ylim = c(0, max(threshold,z))
    plot(z, ylim = ylim, col = "blue2", pch = 16, cex = 1, ylab = "Filter values",
     type = "l", cex.axis = .8, cex.lab = .8)

    points(z, col = "blue2", cex = 0.5, pch = 16)
  }

  filter.loc = apply.filter(fit, filter.bandwidth, threshold)
  rug(filter.loc, col = 1, lwd = 2, ticksize = 0.07,  side = 1)

  rug(enumerate.jumps(fit), col = 1, lwd = 2, ticksize = 0.07, side = 3)

  lines(x = c(-n, 2*n), y = rep(threshold, 2), lty = 2, lwd = 2, col = 2)

  invisible()
  
}
