fem.estep <- function(prms,Y,U){
# require('MASS')
	# Initialization
	Y = as.matrix(Y)
	n = nrow(Y)
	p = ncol(Y)
	K = prms$K
	mu = prms$mean
	prop = prms$prop
	D = prms$D
	b = prms$b
	d = ncol(U)
	QQ = matrix(NA,n,K)
	T = matrix(NA,n,K)
	
	# Compute posterior probabilities
	for (k in 1:K){
	  bk = b[k]
	  mY = prms$my[k,]
	  YY = Y-t(matrix(rep(mY,n),p,n)) 
	  projYY = YY %*% U %*% t(U)
	  if (d==1){
	    QQ[,k] =  1/D[k,1,1] * rowSums(projYY^2) + 1/bk*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k]) + p*log(2*pi)
	  } else{
	    tmp = eigen(D[k,(1:d),(1:d)])
	    A = projYY %*% U %*% tmp$vect %*% diag(sqrt(1/tmp$val))
	    QQ[,k] = rowSums(A^2) + 1/bk*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k]) + p*log(2*pi)
	  }
	}
	# Compute the log-likelihood
	A = -1/2 * QQ
 	loglik = sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
 	
 	# Compute posterior probabilities
	for (k in 1:K) {T[,k] = 1 / rowSums(exp((QQ[,k]*matrix(1,n,K)-QQ)/2))}
	
	# Return the results
	list(T=T,loglik=loglik)
}

