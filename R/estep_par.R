estep_par <-
function(prms,Y,U,memlims,nbcore){
	# Initialization
	Y = as.matrix(Y)
	n = nrow(Y)
	p = ncol(Y)
	K = prms$K
	mu = prms$mean
	prop = prms$prop
	D = prms$D
	d = min(p-1,(K-1))
	QQ = matrix(NA,n,K)
	T = matrix(NA,n,K)
	
	# Compute posterior probabilities (using multicore parallel computing)
	library(doMC)
	registerDoMC(nbcore)
	for (k in 1:K){
	  if (n%%memlim==0) quot = n %/% memlim # quotient de la division euclidienne
	  else quot = n %/% memlim + 1
	  bk = D[k,p,p]
	  mY = prms$my[k,]
	  QQ[,k] = foreach(i=1:quot,.combine='c') %dopar% {
      ni = min(memlim*i,n) - (memlim*(i-1))
	    YY = Y[(memlim*(i-1)+1):min(memlim*i,n),] - t(matrix(rep(mY,ni),p,ni)) 
	    projYY = YY %*% U %*% t(U)
	    if (d==1){
	      QQik =  1/D[k,1,1] * rowSums(projYY^2) + 1/D[k,p,p]*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k])
	    } else{
	      sY = U %*% ginv(D[k,(1:d),(1:d)]) %*% t(U)
	      QQik =   diag(projYY %*% sY %*% t(projYY)) + 1/bk*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k])
	    }
      QQik
	  }
	}
  
	# Compute the log-likelihood
	A = -1/2 * QQ
 	loglik = sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
	for (k in 1:K) {T[,k] = 1 / rowSums(exp(0.5*(QQ[,k]*matrix(1,n,K)-QQ)))}
	
	# Return the results
	list(T=T,loglik=loglik)
}

