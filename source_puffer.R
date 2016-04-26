#using the puffer transformation on the fused lasso
puffer.estimate <- function(y){
  n = length(y)

  #form Y  
  y.centered = y - mean(y)
 
  #form X
  X = matrix(0, nrow = n, ncol = n-1)
  X[lower.tri(X, diag = F)] = 1
  X.centered = scale(X, center = T, scale = F)

  #form A
  A = matrix(0, n, n)
  A[lower.tri(A, diag = T)] = 1
  A.inv = solve(A)

  #form F
  svd.decomp = svd(X.centered)
  F.mat = svd.decomp$u %*% diag(1/svd.decomp$d) %*% t(svd.decomp$u)

  #transform X.centered into Z
  Z = F.mat %*% X.centered
 
  #transform y.centered into a
  a = F.mat %*% y.centered

  #now run lasso
  res = cv.glmnet(x = Z, y = a, intercept = F)
  res.coef = as.numeric(coef(res$glmnet.fit, s = res$lambda.1se)[-1])

  #transform back to fused lasso solution
  first.val = mean(a) - apply(Z, 2, mean)%*%res.coef
  res.coef = c(first.val, res.coef)

  #WARNING: not sure if I need to add mean(y) back in??
  A%*%res.coef
}
