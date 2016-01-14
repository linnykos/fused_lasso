
#based on the recursive formulation in ryan's paper, output the numerator
#  and denominator for each jump
recursive.solve <- function(y,num.est.jumps=5){
  n = length(y)
  res.list = list()
  lambda.vec = numeric(num.est.jumps)
  D = form.diff.matrix(n)
  total.index = 1:n
  
  #do the first iteration manually
  Z = solve(D%*%t(D))%*%D
  num = Z%*%y
  lambda.vec[1] = max(abs(num))
  idx = which.max(abs(num))
  s = sign(num[idx])
  total.index = total.index[-idx]
  res.list[[1]] = list(s = s, idx = idx, lambda = lambda.vec[1], num = as.vector(num), denom = sign(num))
  
  for(iter in 2:num.est.jumps){
    Z2 = solve(D[-idx,]%*%t(D[-idx,]))%*%D[-idx,]
    num = Z2%*%y
    denom = Z2%*%t(D[idx,,drop=FALSE])%*%s
    tmp1 = denom+1
    tmp2 = denom-1
    tmp1[which(tmp1<1e-5)] = NA
    tmp2[which(tmp1<1e-5)] = NA
    val = apply(cbind(num/tmp1,num/tmp2),1,max,na.rm=TRUE)
    lambda.vec[iter] = max(val)
    idx = c(idx, total.index[which.max(val)])
    stmp = if(tmp[which.max(val)]+1>0) 1 else -1
    s = c(s,stmp)
    total.total = total.index[-which.max(val)]
    
    
    #reorder
    tmporder = order(idx)
    idx = idx[tmporder]
    s = s[tmporder]
    res.list[[iter]] = list(s = s, idx = idx, lambda = lambda.vec[iter], 
                            num = as.vector(num), denom = as.vector(denom))
  }
  
  res.list
}