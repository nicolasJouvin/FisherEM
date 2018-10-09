fem.criteria <-
function(loglik,T,prms,n){
	K = prms$K
	d = K-1
	V = prms$V
	if (prms$method=='sparse'){ p = sum(abs(V) > 1e-2) } else {p = prms$p}
	comp = switch(as.character(prms$model),
		'DkBk' = (K-1) + K*d + (K-1)*(p-K/2) + K^2*(K-1)/2 + K,
		'DkB'  = (K-1) + K*d + (K-1)*(p-K/2) + K^2*(K-1)/2 + 1,
		'DBk'  = (K-1) + K*d + (K-1)*(p-K/2) + K*(K-1)/2 + K,
		'DB'   = (K-1) + K*d + (K-1)*(p-K/2) + K*(K-1)/2 + 1,
		'AkjBk'= (K-1) + K*d + (K-1)*(p-K/2) + K^2,
		'AkjB' = (K-1) + K*d + (K-1)*(p-K/2) + K*(K-1)+1,
		'AkBk' = (K-1) + K*d + (K-1)*(p-K/2) + 2*K,
		'AkB'  = (K-1) + K*d + (K-1)*(p-K/2) + K+1,
		'AjBk' = (K-1) + K*d + (K-1)*(p-K/2) + (K-1)+K,
		'AjB'  = (K-1) + K*d + (K-1)*(p-K/2) + (K-1)+1,
		'ABk'  = (K-1) + K*d + (K-1)*(p-K/2) + K+1,
		'AB'   = (K-1) + K*d + (K-1)*(p-K/2) + 2)
	aic = loglik - comp # AIC criterion
	bic = loglik - 1/2 * comp * log(n) # BIC criterion
	Z = ((T - apply(T,1,max))==0) + 0
	icl = loglik - 1/2 *  comp * log(n) - sum(Z*log(T+1e-15)) # ICL criterion
	list(aic=aic,bic=bic,icl=icl,comp=comp)
}

