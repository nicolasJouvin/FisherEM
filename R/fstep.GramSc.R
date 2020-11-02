fstep.GramSc <-
function(XX,d,T,S,kernel){
	n = nrow(XX)
	p = ncol(XX)
	K = ncol(T)
	#m = colMeans(Y)
	if (is.null(d))	d = min(p-1,(K-1))
	U = matrix(NA,p,d)

	# Compute S
	#XX = as.matrix(Y - t(m*t(matrix(1,n,p))))
	TT = t(apply(T,1,"/",sqrt(colSums(T))))
	
	if (n>p & kernel==''){
		#S = cov(Y)*(n-1)/n
		B = t(XX)%*%TT%*%t(TT)%*%XX / n
	
		# Eigendecomposition of S^-1 * B
		eig = eigen(ginv(S)%*%B)
		#eig = eigs(ginv(S)%*%B,1) # requires the rARPACK library
		U[,1] = matrix(Re(eig$vec[,1]),ncol=1)

		if (d>1){
		  for (k in 2:d){ 
		    W = as.matrix(U[,-c(k:d)])  
		    
		    # start of the gram-schmidt algorithm
		    base = diag(p)
		    base <- as.matrix(base[,1:(p-k+1)])
		    W =  cbind(W,base)
		    v <- W
		    for (l in 2:p){
		      proj <- c(crossprod(W[,l],v[,1:(l-1)])) / diag(crossprod(v[,1:(l-1)]))
		      v[,l] <- matrix(W[,l],ncol=1) - matrix(v[,1:(l-1)],nrow=p) %*%  matrix(proj,ncol=1)
		      v[,l] <- v[,l] /  rep(sqrt(crossprod(v[,l])),p)
		    }
		    P = v[,k:ncol(v)]
		    
		    # P = qr.Q(qr(W),complete=TRUE)[,k:p]
		    B.p = crossprod(P,crossprod(t(B),P))
		    S.p = crossprod(P,crossprod(t(S),P))
		    eig = eigen(ginv(S.p)%*%B.p)
		    #eig = eigs(ginv(S.p)%*%B.p,1) # requires the rARPACK library
		    
		    Umax = matrix(Re(eig$vec[,1]),ncol=1) 
		    U[,k]= P%*%Umax	
		  }
		}
	}
	if (d==1) {U = U[,1]} 
	as.matrix(U)
}

