#more of just simple experiments than tests
rm(list=ls())
source("source_header.R")

sigma = 1
n.length = 10
n.vec = round(10^seq(2, 4, length.out = n.length))
trials = 50
jump.mean =     c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)

i = 6
trial = 1
truth = form.truth(jump.mean, jump.location, n.vec[i])
true.jumps = enumerate.jumps(truth)

set.seed(i*trial*10)
    
y = generate.problem(n.vec[i], jump.mean, jump.location, sigma)

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

library(glmnet)
set.seed(10)
res = cv.glmnet(x = X, y = y)
res.coef = as.numeric(coef(res$glmnet.fit, s = res$lambda.1se))
fit = A%*%res.coef

set.seed(10)
res2 = cv.glmnet(x = X.centered, y = y.centered, intercept = F)
res.coef2 = as.numeric(coef(res2$glmnet.fit, s = res2$lambda.1se)[-1])
first.val = mean(y) - apply(X, 2, mean)%*%res.coef2
res.coef2 = c(first.val, res.coef2)
fit2 = A%*%res.coef2

assert_that(abs(sum(fit - fit2)) < 1e-6)

#form F
svd.decomp = svd(X.centered)
F.mat = svd.decomp$u %*% diag(1/svd.decomp$d) %*% t(svd.decomp$u)

#transform X.centered into Z
Z = F.mat %*% X.centered

#transform y.centered into a
a = F.mat %*% y.centered

tmp = svd(Z)
assert_that(all(tmp$d == 1))

res3 = cv.glmnet(x = Z, y = a, intercept = F)
res.coef3 = as.numeric(coef(res3$glmnet.fit, s = res3$lambda.1se)[-1])
first.val = mean(a) - apply(Z, 2, mean)%*%res.coef3
res.coef3 = c(first.val, res.coef3)
fit3 = A%*%res.coef3


