convex_clust <- function(y, lambda1, lambda2,
	rho = 1, init.beta = NA, init.alpha = NA, init.u = NA,
	tol = 0.01, max.iter = 25, verbose = FALSE){
	#do ADMM
	n = length(y)
	if(missing(init.beta)) init.beta = y
	if(missing(init.alpha)) init.alpha = y
	if(missing(init.u)) init.u = rep(0,n)

	obj.val = 10^8
	prev.obj.val = 10^7
	beta = init.beta
	alpha = init.alpha
	u = init.u
	iter = 1

	while(abs(obj.val-prev.obj.val)>tol ){
	  prev.obj.val = obj.val
    
    if(verbose){
      obj.val = .25*sum((y-beta)^2) + lambda1*sum(abs(diff(beta)))+.25*sum((y-alpha)^2)+lambda2*sum(sapply(2:n,function(x){sum(abs(alpha[x]-alpha[1:(x-1)]))}))
      print(paste("Iteration ", iter," start: ",obj.val))
      
    }
	  
		#opt over beta, trend filer
		y.tmp = (y+(rho/2)*(u+alpha))/(1+rho/2)
		out = fusedlasso1d(y.tmp)
		beta = as.numeric(coef.genlasso(out, lambda = lambda1*2)$beta)
    
	  if(verbose){
	    obj.val = .25*sum((y-beta)^2) + lambda1*sum(abs(diff(beta)))+.25*sum((y-alpha)^2)+lambda2*sum(sapply(2:n,function(x){sum(abs(alpha[x]-alpha[1:(x-1)]))}))
	    print(paste("beta: ",obj.val))
	  }

		#opt over alpha, convex clus
		y.tmp = (y+(rho/2)*(beta-u))/(1+rho/2)
		y.tmp = t(y.tmp)
	  w <- kernel_weights(y.tmp,0.5)
	  w <- knn_weights(w,5,n)
		out = suppressWarnings(cvxclust(y.tmp, w = w, gamma = lambda2*2))
		alpha = as.numeric(out$U[[1]])
    alpha = round(alpha,2)
    
	  if(verbose){
	    obj.val = .25*sum((y-beta)^2) + lambda1*sum(abs(diff(beta)))+.25*sum((y-alpha)^2)+lambda2*sum(sapply(2:n,function(x){sum(abs(alpha[x]-alpha[1:(x-1)]))}))
	    print(paste("alpha: ",obj.val))
      print(paste("K = ",length(table(alpha))))
	  }
    
		#opt over u
		u = u + alpha - beta
	  if(verbose){
      print(paste("u: ",sum(abs(alpha-beta))))
	  }
	
		
		obj.val = .25*sum((y-beta)^2) + lambda1*sum(abs(diff(beta)))+.25*sum((y-alpha)^2)+lambda2*sum(sapply(2:n,function(x){sum(abs(alpha[x]-alpha[1:(x-1)]))}))
	
		if(iter > max.iter ) break()
		iter = iter +1
    if(iter%%floor(max.iter/10)==0) cat('*')
    
    if(iter > 5 & sum(abs(alpha-beta))/n < tol*4) break()

	}

  K = length(table(alpha))
  
	list(beta = beta, alpha = alpha, iter = iter, obj.val = obj.val, u = u, K = K)
}
