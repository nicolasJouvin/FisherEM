fstep.fisher <- function(XX,d,T,S,kernel){
	n = nrow(XX)
	p = ncol(XX)
	K = ncol(T)
	#m = colMeans(Y)
	if (is.null(d))	d = min(p-1,(K-1))
	
	# Compute S
	#XX = as.matrix(Y - t(m*t(matrix(1,n,p))))
	TT = t(apply(T,1,"/",sqrt(colSums(T))))
	
	if (n>p & kernel==''){
		#S = t(XX) %*% XX /n
		B = t(XX)%*%TT%*%t(TT)%*%XX / n
	
		# Eigendecomposition of S^-1 * B
		eig = svd(ginv(S)%*%B,nu=d,nv=0)
		#eig = svds(ginv(S)%*%B,d,nu=d,nv=0)
		U   = eig$u[,1:d]
	}
	else{
		if (n<p | kernel=='linear') G = XX %*% t(XX)
		if (kernel=='rbf') {sigma=1; G = as.matrix(exp(dist(XX,diag=T)^2/(2*sigma^2)))}
		if (kernel=='sigmoid') {a=1;r=0.1;G = tanh(a * XX %*% t(XX) + r)}
		lambda = 0.5
		S = G %*% G + lambda*diag(n)
		B = G %*% TT %*% t(TT) %*% G
		H = svd(ginv(S)%*%B,nu=d,nv=0)$u[,1:d]
		U = svd(t(XX) %*% H,nu=d,nv=0)$u[,1:d]
	}
  return(as.matrix(U))
}

