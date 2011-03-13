## Fisher EM ##
#
# Author: 	C.Bouveyron & C.Brunet
# Institute:	SAMM - Paris 1 / Evry
# Date: 	February 2011
# License:	GPL v2
# 
# ***************************************************************
# Main function which computes a discriminant subspace 
# alternatively with Gaussian Mixture Model
# ***************************************************************
# Call of libraries needed
library('MASS')
# ***************************************************************
fem <- function(Y,K,init='random',maxit=100,eps=1e-6,Tinit=c(),model='AkjBk',kernel=''){

    Mod <- c("AkjBk","AkjB", "AkBk", "AkB", "AjBk", "AjB", "ABk", "AB", "DkBk", "DkB",
             "DBk", "DB")
    Init <- c("random", "kmeans", "mini-em", "param")
    if (any(Mod == model) == 0) 
        stop("Error : invalid model name\n", .Call = FALSE)
    else if (!is.numeric(K) || min(K) < 1) 
        cat("Error : K must be a positive integers\n")
    else if (K>=ncol(Y)) 
        cat("Error : the number of K is incorrect\n")

		
# Usage: res <- fem(Y,K,...){
	# Initialization
	Y = as.matrix(Y)
	n = nrow(Y)
	p = ncol(Y)
	d = min((K-1),(p-1))
# 
	# New objects
	Lobs = rep(c(-Inf),1,(maxit+1))
# 	
	# Initialization of T
	if (init=='user'){ T = Tinit}
	if (init=='kmeans'){
		T = matrix(0,n,K)
		ind = kmeans(Y,K)$cluster
		for (i in 1:n){ T[i,ind[i]] = 1 }
	}
	if (init=='random'){
		T = t(rmultinom(n,1,c(rep(1/K,K))))
	}
	if (init=='param'){
		V = princomp(Y)$loadings[1:d,]
		X = as.data.frame(as.matrix(Y) %*% V[,1:d])
		prop = rep(c(1/K),K)
		mu = mvrnorm(K,colMeans(X),cov(X))
		my = mu %*% t(V) # selon le modele 	
		D = matrix(NA,K,p)
		res = eigen(cov(X)) 
		D[,1:d] = mean(res$val) # il n y a que d valeurs propres 
		D[,(d+1):p] = mean(eigen(cov(Y))$val[(d+1):p]) 	
		prms = list(K=K,p=p,mean=mu,my=my,prop=prop,V=V,D=D)
		T = estep(prms,Y,prms$V)
	}
	if (init=='miniEM'){
		source('../EM.r')
		TT = array(NA,c(10,n,K))
		Lik = rep(NA,10)
		for (nb in 1:10){
			res.em = EM(Y,K,model='full',maxit=10)
			Lik[nb] = res.em$lik
			TT[nb,,] = res.em$P
		}
		T = TT[which.max(Lik),,]
	}

	V         = fstep(Y,T,kernel)
	prms      = mstep(Y,V,T,model=model)
	res.estep = estep(prms,Y,V)
	T         = res.estep$T
	Lobs[1]   = res.estep$loglik

	# Main loop
	Linf_new  = Lobs[1]
	for (i in 1:maxit){
		cat('*')

		# The three main steps F, M and E
		V         = fstep(Y,T,kernel)
		prms      = mstep(Y,V,T,model=model)
		res.estep = estep(prms,Y,V)
		T         = res.estep$T
		Lobs[i+1] = res.estep$loglik
		
		# Stop criterion
		if (i>=2){
			acc      = (Lobs[i+1] - Lobs[i]) / (Lobs[i] - Lobs[i-1])
			Linf_old = Linf_new
			Linf_new = Lobs[i] + 1/(1-acc) * (Lobs[i+1] - Lobs[i])
			if (abs(Linf_new - Linf_old) < eps) {break}
		}

	}
	cat('\n')
	
	# Returning the results
	cls  = max.col(T)
	crit = criteria(Lobs[(i+1)],T,prms,n); # BIC
	res  = list(cls=cls,P=T,prms=prms,V=V,aic=crit$aic,bic=crit$bic,icl=crit$icl,loglik=Lobs[2:(i+1)],ll=Lobs[i+1],data=Y)
	class(res)="fem"
	res
	
}


########################################################################
# FISHER STEP:  computes the orthonormal disciminant subspace owing to the Fisher criterion.
# Method used : Singular value decomposition
#-------------------------------------------------------------------
fstep <- function(Y,T,kernel){
	n = nrow(Y)
	p = ncol(Y)
	K = ncol(T)
	m = colMeans(Y)
	d = min(p-1,(K-1))

	# Compute S
# browser()
	XX = as.matrix(Y - t(m*t(matrix(1,n,p))))
	TT = t(apply(T,1,"/",sqrt(colSums(T))))
	
	if (n>p & kernel==''){
		S = t(XX) %*% XX /n
		B = t(XX)%*%TT%*%t(TT)%*%XX / n
		
		# Eigendecomposition of S^-1 * B
		eig = svd(ginv(S)%*%B,nu=d,nv=0)
		U = eig$u[,1:d]
	}
	else{
		cat('Kernel mode!\n')
		if (n<p | kernel=='linear') G = XX %*% t(XX)
		if (kernel=='rbf') {sigma=1; G = as.matrix(exp(dist(XX,diag=T)^2/(2*sigma^2)))}
		if (kernel=='sigmoid') {a=1;r=0.1;G = tanh(a * X %*% t(X) + r)}
		lambda = 0
		S = G %*% G + lambda*diag(n)
		B = G %*% TT %*% t(TT) %*% G
		H = svd(ginv(S)%*%B,nu=d,nv=0)$u[,1:d]
		V = svd(t(Y) %*% H,nu=d,nv=0)$u[,1:d]
	}
}
############################################################################
# E-STEP : computes the posterior probabilities of each observations
#-------------------------------------------------------------------
estep <- function(prms,Y,U){
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
	memlim = 2500
	
	# Compute posterior probabilities
	for (k in 1:K){
		bk = D[k,p,p]
		mY = prms$my[k,]
		YY = Y-t(matrix(rep(mY,n),p,n)) 
		projYY = YY %*% U %*% t(U)
		if (nrow(YY)<=memlim){
			if (d==1){
				QQ[,k] =  1/D[k,1,1] * rowSums(projYY^2) + 1/D[k,p,p]*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k])
			} else{
				sY = U %*% ginv(D[k,(1:d),(1:d)]) %*% t(U)
				QQ[,k] =   diag(projYY %*% sY %*% t(projYY)) + 1/bk*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k])
			}
		}
		else{ # Do the same but by blocks of memlim obs.
			quot = n %/% memlim # quotient de la division euclidienne
			for (i in 1:(quot+1)){
				if (d==1){
					QQ[(memlim*(i-1)+1):min(memlim*i,n),k] =  1/D[k,1,1] * rowSums(projYY[(memlim*(i-1)+1):min(memlim*i,n),]^2) + 1/D[k,p,p]*rowSums((YY[(memlim*(i-1)+1):min(memlim*i,n),] - projYY[(memlim*(i-1)+1):min(memlim*i,n),])^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k])
				} else{
					sY = U %*% ginv(D[k,(1:d),(1:d)]) %*% t(U)
					QQ[(memlim*(i-1)+1):min(memlim*i,n),k] =   diag(projYY[(memlim*(i-1)+1):min(memlim*i,n),] %*% sY %*% t(projYY[(memlim*(i-1)+1):min(memlim*i,n),])) + 1/bk*rowSums((YY[(memlim*(i-1)+1):min(memlim*i,n),] - projYY[(memlim*(i-1)+1):min(memlim*i,n),])^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k])
				}
			}
		
		}
	}
	# Compute the log-likelihood
	A = -1/2 * QQ
 	loglik = sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
	for (k in 1:K) {T[,k] = 1 / rowSums(exp(0.5*(QQ[,k]*matrix(1,n,K)-QQ)))}
	
	# Return the results
	list(T=T,loglik=loglik)
}
#######################################################################################
# M-STEP : updates the parameters of the Mixture model (mean, prior probabilities ...)
#-------------------------------------------------------------------------------------
mstep <- function(Y,U,T,model){
	# 12 different submodels: [DkBk] ... [AkjBk]
	# Initialization
	Y = as.matrix(Y)
	n = nrow(Y)
	p = ncol(Y)
	K = ncol(T)
	d = min(p-1,(K-1))

	mu   = matrix(NA,K,K-1)
	m   = matrix(NA,K,p)
	prop = rep(c(NA),1,K)
	D = array(0,c(K,p,p))

	# Projection
	X = Y %*% U

	# Estimation
	for (k in 1:K){

		nk  = sum(T[,k])
		# Prior Probability
		prop[k] = nk / n
		# Mean in the latent space
		mu[k,]  = colSums((T[,k]*matrix(1,n,d))* X) / nk
		# Observed space
		m[k,]  = colSums(T[,k]*matrix(1,n,p)* Y) / nk
		YY  = as.matrix(Y - t(m[k,]*matrix(1,p,n)))
		Ck  = crossprod(T[,k]*matrix(1,n,p)* YY, YY) / (nk-1) #crossprod(x,y) = t(x) %*% y
		C   = cov(Y)

		# Estimation of Delta k amongst 8 submodels
		if (model=='DkBk'){
			D[k,(1:d),(1:d)] = crossprod(Ck%*%U,U)
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		if (model=='DkB'){
			D[k,(1:d),(1:d)] = crossprod(Ck%*%U,U)
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		if (model=='DBk'){
			D[k,(1:d),(1:d)] = crossprod(C%*%U,U)
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		if (model=='DB'){
			D[k,(1:d),(1:d)] = crossprod(C%*%U,U)
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		
		if (model=='AkjBk'){
		if (d==1){D[k,1,1] = diag(crossprod(Ck%*%U,U))} else {
			D[k,(1:d),(1:d)] = diag(diag(crossprod(Ck%*%U,U)))}
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		
		if (model=='AkjB'){
		if (d==1){D[k,1,1] = diag(crossprod(Ck%*%U,U))} else {
			D[k,(1:d),(1:d)] = diag(diag(crossprod(Ck%*%U,U)))}
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		
		if (model=='AkBk'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(Ck%*%U,U)))/d} else{
			D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(Ck%*%U,U)))/d,d))}
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		
		if (model=='AkB'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(Ck%*%U,U)))/d} else{
			  D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(Ck%*%U,U)))/d,d))}
			  bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			  bk[bk<=0] = 1e-3
			  D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		
		if (model=='AjBk'){
		if (d==1){D[k,1,1] = diag(crossprod(C%*%U,U))} else {
			D[k,(1:d),(1:d)] = diag(diag(crossprod(C%*%U,U)))}
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		if (model=='AjB'){
		if (d==1){D[k,1,1] = diag(crossprod(C%*%U,U))} else{
			D[k,(1:d),(1:d)] = diag(diag(crossprod(C%*%U,U)))}
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		if (model=='ABk'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(C%*%U,U)))} else {
			D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(C%*%U,U)))/d,d))}
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
		}
		if (model=='AB'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(C%*%U,U)))} else {
			D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(C%*%U,U)))/d,d))}
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))	
		}
	}
	prms = list(K=K,p=p,mean=mu,my=m,prop=prop,D=D,model=model)
}
######################################################################################
# Function which computes 3 different criteria to choose the
# number of classes:
#		* AIC : Akaike Information criterion
# 		* BIC : Bayesian Information criterion
#		* ICL : 
#-------------------------------------------------------------------------------------
criteria <- function(loglik,T,prms,n){
	K = prms$K
	p = prms$p
	comp = switch(prms$model,
		'DkBk' = (K-1) + K*p+ (K-1)*(p-K/2) + K^2*(K-1)/2 + K,
		'DkB'  = (K-1) + K*p+ (K-1)*(p-K/2) + K^2*(K-1)/2 + 1,
		'DBk'  = (K-1) + K*p+ (K-1)*(p-K/2) + K*(K-1)/2 + K,
		'DB'   = (K-1) + K*p+ (K-1)*(p-K/2) + K*(K-1)/2 + 1,
		'AkjBk'= (K-1) + K*p + (K-1)*(p-K/2) + K^2,
		'AkjB' = (K-1) + K*p + (K-1)*(p-K/2) + K*(K-1)+1,
		'AkBk' = (K-1) + K*p + (K-1)*(p-K/2) + 2*K,
		'AkB'  = (K-1) + K*p + (K-1)*(p-K/2) + K+1,
		'AjBk' = (K-1) + K*p + (K-1)*(p-K/2) + (K-1)+K,
		'AjB'  = (K-1) + K*p + (K-1)*(p-K/2) + (K-1)+1,
		'ABk'  = (K-1) + K*p + (K-1)*(p-K/2) + K+1,
		'AB'   = (K-1) + K*p + (K-1)*(p-K/2) + 2)
	aic = loglik - comp # AIC criterion
	bic = loglik - 1/2 * comp * log(n) # BIC criterion
	icl = loglik - 1/2 *  comp * log(n) - sum(T*log(T)) # ICL criterion
	list(aic=aic,bic=bic,icl=icl)
}