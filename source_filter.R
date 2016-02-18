apply.filter <- function(fit, bandwidth, threshold, y = NA, return.type = c("location", "filter", "refit")){
  n = length(fit)
  return.type = return.type[1]
  
  if(return.type == "refit") assert_that(all(!is.na(y)))
  
  jump.loc = enumerate.jumps(fit)

  z = sapply((bandwidth+1):(n-bandwidth+1), function(x){
    mean(fit[x:(x+bandwidth-1)]) - mean(fit[(x-bandwidth):(x-1)])
  })
  
  #pad the z's
  z = c(rep(0, bandwidth), z, rep(0, bandwidth-1))
  if(return.type == "filter") return(z)
  
  z.thresIdx = which(abs(z) >= threshold)
  
  jump.remain = intersect(jump.loc, z.thresIdx)
  if(return.type == "location") return(jump.remain)
  
  jump.remain = c(1, jump.remain, n+1)
  refit = sapply(2:length(jump.remain), function(x){
    mean(y[x-1]:(y[x]-1))
  })
  
  return(refit)
}