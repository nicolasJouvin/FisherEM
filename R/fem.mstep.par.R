fem.mstep.par <- function(Y,U,T,model,method){
	# 12 different submodels: [DkBk] ... [AkjBk]
	# Initialization
	Y = as.matrix(Y)
	n = nrow(Y)
	p = ncol(Y)
	K = ncol(T)
	d = min(p-1,(K-1))
	U = as.matrix(U)

	mu   = matrix(NA,K,d)
	m   = matrix(NA,K,p)
	prop = rep(c(NA),1,K)
	D = array(0,c(K,d,d))
	b = rep(NA,K)

	# Projection
	X = as.matrix(Y) %*% as.matrix(U)

	# Estimation
	test = 0
	fun <- function(k){
		nk  = sum(T[,k])
		if (nk ==0) stop("some classes become empty\n",call.=FALSE)
		# Prior Probability
		prop[k] = nk / n
		# Mean in the latent space
		mu[k,]  = colSums((T[,k]*matrix(1,n,d))* X) / nk
		# Observed space
		m[k,]  = colSums(T[,k]*matrix(1,n,p)* Y) / nk
		YY  = as.matrix(Y - t(m[k,]*matrix(1,p,n)))
		if (nk < 1) denk = 1 else denk = (nk-1)
		YYk = T[,k] * matrix(1,n,p) * YY
		#browser()
		if (model %in% c('DkB','DB','AkjB','AkB','AjB','AB')){
		  TT = t(apply(T,1,"/",sqrt(colSums(T))))
		  W = (crossprod(YY) - crossprod(t(TT)%*%YY)) / n
		}
		
		# Estimation of Delta k amongst 8 submodels
		if (model=='DkBk'){ # OK !!!
			D[k,(1:d),(1:d)] = crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk
			bk = (sum(YYk*YY)/denk - sum(diag(D[k,(1:d),(1:d)]))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='DkB'){
			D[k,(1:d),(1:d)] =  crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk
			bk = (sum(diag(W)) - sum(diag(crossprod(W%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='DBk'){
			D[k,(1:d),(1:d)] = crossprod(YY%*%U,YY%*%U)/n
			bk = (sum(YYk*YY)/denk - sum(diag(D[k,(1:d),(1:d)]))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='DB'){
			D[k,(1:d),(1:d)] = crossprod(YY%*%U,YY%*%U)/n
			bk = (sum(YY*YY)/n - sum(diag(D[k,(1:d),(1:d)]))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='AkjBk'){
		if (d==1){D[k,1,1] = diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk)} else {
			D[k,(1:d),(1:d)] = diag(diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk))}
			bk = (sum(YYk*YY)/denk - sum(diag(D[k,(1:d),(1:d)]))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='AkjB'){
		if (d==1){D[k,1,1] = diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk)} else {
			D[k,(1:d),(1:d)] = diag(diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk))}
			bk = (sum(YY*YY)/n - sum(diag(D[k,(1:d),(1:d)]))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='AkBk'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk))/d} else{
			D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk))/d,d))}
			bk = (sum(YYk*YY)/denk - sum(diag(D[k,(1:d),(1:d)]))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='AkB'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk))/d} else{
			  D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk))/d,d))}
			  bk = (sum(YY*YY)/n - sum(diag(D[k,(1:d),(1:d)]))) / (p-d)
			  bk[bk<=0] = 1e-3
			  b[k] = bk
   			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='AjBk'){
		if (d==1){D[k,1,1] = diag(crossprod(YY%*%U,YY%*%U)/n)} else {
			D[k,(1:d),(1:d)] = diag(diag(crossprod(YY%*%U,YY%*%U)/n))}
			bk = (sum(YYk*YY)/denk - sum(diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='AjB'){
		if (d==1){D[k,1,1] = diag(crossprod(YY%*%U,YY%*%U)/n)} else{
			D[k,(1:d),(1:d)] = diag(diag(crossprod(YY%*%U,YY%*%U)/n))}
			bk = (sum(YY*YY)/n  - sum(diag(D[k,(1:d),(1:d)]))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='ABk'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(YY%*%U,YY%*%U)/n))} else {
			D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(YY%*%U,YY%*%U)/n))/d,d))}
			bk = (sum(YYk*YY)/denk - sum(diag(crossprod(T[,k]*matrix(1,n,d)* YY %*%U, YY%*%U) / denk))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		else if (model=='AB'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(YY%*%U,YY%*%U)/n))} else {
			D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(YY%*%U,YY%*%U)/n))/d,d))}
			bk = (sum(YY*YY)/n  - sum(diag(D[k,(1:d),(1:d)]))) / (p-d)
			bk[bk<=0] = 1e-3
			b[k] = bk
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
	}
	prms = list(K=K,p=p,mean=mu,my=m,prop=prop,D=D,b=b,model=model,method=method,V=U,test=test)
}

