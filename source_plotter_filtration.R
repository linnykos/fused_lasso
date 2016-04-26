plot.filter.filtration <- function(fit, filter.bandwidth, return.plot = T,
 return.val = F, discretization = 100){
  filter.val = apply.filter(fit, filter.bandwidth, return.type = "filter")

  n = length(fit)

  max.val = max(filter.val)
  min.val = 0
  
  thres.val = seq(min.val, max.val, length.out = discretization)
  idx.list = vector("list", length(thres.val))

  if(return.plot) plot(NA, xlim = c(1, n), ylim = c(min.val, max.val),
   main = "Screened Changepoints Filtration", xlab = "Index", ylab = "Filter Threshold")

  for(i in 1:length(thres.val)){
    idx.list[[i]] = apply.filter(fit, filter.bandwidth, threshold = thres.val[i], 
     return.type = "location")

    if(return.plot & length(idx[[i]]) > 0) sapply(idx.list[[i]], function(x){
       points(x = x, y = thres.val[i], pch = 16)
       invisible()
     })
  }

  if(return.val) idx.list else invisible()
}

#multiscale version of the filter filtration
plot.filter.filtratrion.multiscale <- function(){
}
