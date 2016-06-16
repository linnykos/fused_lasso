
.plot.helper <- function(jump.location, jump.mean, n, col = "black", 
 lwd = 3, lty = 1, plot.vertical = T){
  assert_that(length(jump.location) == length(jump.mean))
  assert_that(all(jump.location == sort(jump.location, decreasing = F)))

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
                      filter.bandwidth = NA, threshold = NA){
  
  n = length(fit)
  true.seq = form.truth(jump.mean, jump.location, n)
  if(is.na(mse)) mse = compute.mse(fit, true.seq = true.seq)
  if(is.na(num.est.jumps)) num.est.jumps = length(enumerate.jumps(fit))
 
  par(mfrow=c(2,1),mar=c(2,4,1,1))
  
  plot(y, col=rgb(.7,.7,.7), pch=16, cex=1.25, ylab = "Data values", 
   cex.axis = .8, cex.lab = .8)

  .plot.primal(jump.mean, jump.location, y, fit, tol)
 
  if(is.na(filter.bandwidth)) filter.bandwidth = ceiling(0.25*log(n)^2)
    
  .plot.filter(fit, filter.bandwidth, jump.mean, jump.location, lambda, mse, 
    threshold, num.est.jumps)
 
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
 fit, tol, truth = NA, verbose = F, truth.col = 2, est.col = 4){
  n = length(y)
 
  if(all(is.na(truth))){
    true.seq = form.truth(jump.mean, jump.location, n)
  
    #plot truth
    tmp = seq(0, 1,length.out = n)
    jump.location2 = .extract.location(jump.location, tmp)
    .plot.helper(jump.location2, jump.mean, n, col = 2)
  } else {
    truth.split = split.signal(truth)
    .plot.helper(truth.split$location, truth.split$mean, n, col = truth.col)
  }

  res.split = split.signal(fit)
  .plot.helper(res.split$location, res.split$mean, n, col = est.col)
  
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

#WARNING: Fix the default for num.est.jumps
.plot.filter <- function(fit, filter.bandwidth, jump.mean, jump.location, 
 lambda, mse, threshold = NA, num.est.jumps = NA, verbose = F){

  n = length(fit)
  
  if(is.na(threshold)) threshold = min(abs(diff(jump.mean)))/2
  
  z = apply.filter(fit, filter.bandwidth, threshold, return.type = "filter")
 
  ylim = c(0, max(threshold,z))

  ylim[1] = ylim[1] - 0.1*diff(range(ylim))
  ylim[2] = ylim[2] + 0.1*diff(range(ylim))

  plot(z, ylim = ylim, col = 4, pch = 16, cex = 1, 
   ylab = "Absolute filter values", type = "l", cex.axis = .8, cex.lab = .8)

  points(z, col = 4, cex = 0.5, pch = 16)
  
  filter.loc = apply.filter(fit, filter.bandwidth, threshold)
  rug(filter.loc, col = 1, lwd = 2, ticksize = 0.07,  side = 3)

  rug(enumerate.jumps(fit), col = 1, lwd = 2, ticksize = 0.07, side = 1)

  lines(x = c(-n, 2*n), y = rep(threshold, 2), lty = 2, lwd = 2, col = 3)

  invisible()
}
