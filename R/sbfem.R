sbfem <- function(Y,K=2:6,obj=NULL,model='AkjBk',method='reg',crit='icl', maxit.em=100,
                 eps.em=1e-6, maxit.ve=3, eps.ve=1e-4, lambda = 1e3, emp.bayes=T,
                 init='kmeans', nstart=5,Tinit=c(),kernel='',disp=FALSE,
                 maxit.sparse=15, eps.sparse=NULL, l1=0.1,l2=0,nbit=2){
  
  # Initial tests & settings
  call = match.call()
  if (length(l2)>1) stop("The l2 parameter must be a scalar value (model selection for the l2 penalty is not yet implemented!).\n",call.=FALSE) 
  if (max(l1)>1 | min(l1)<0 | length(l2)>1 | l2>1 | l2<0) stop("Parameters l1 and l2 must be within [0,1]\n",call.=FALSE)
  bic = aic = icl = c()
  RES = list()
  # Initialization with FEM
  if (is.null(obj)){
    try(res0 <- bfem(Y, K=K, model=model, method=method, crit='icl', 
                     maxit.em=maxit.em, eps.em = eps.em, 
                     maxit.ve=maxit.ve, eps.ve=eps.ve, 
                     lambda = lambda, emp.bayes=emp.bayes,
                     init=init, nstart=nstart, Tinit=Tinit, kernel=kernel, 
                     disp=disp, mc.cores = 1))
  } else{
    if (!is(obj, 'bfem')) stop('An object of class "bfem" is required!')
    else res0 = obj
  }
  
  # Sparse FEM
  if(is.null(eps.sparse)) eps.sparse = eps.em
  
  for (i in 1:length(l1)){
    try(RES[[i]] <- bfem.sparse(Y,res0$K,model=res0$model,
                                maxit.sparse=maxit.sparse,
                                eps.sparse=eps.sparse,
                                maxit.ve=maxit.ve,
                                eps.ve=eps.ve,
                               Tinit=res0$P, 
                               lambda=res0$hypers$lambda,
                               # lambda = lambda,
                               nu = res0$hypers$nu,
                               emp.bayes = emp.bayes,
                               l1=l1[i],l2=l2,nbit=nbit))
    try(bic[i] <- RES[[i]]$bic)
    try(aic[i] <- RES[[i]]$aic)
    try(icl[i] <- RES[[i]]$icl)
  }
  
  # Post-treatment and returning results
  if (crit=='bic'){ id_max = which.max(bic); crit_max = RES[[id_max]]$bic}
  if (crit=='aic'){ id_max = which.max(aic); crit_max = RES[[id_max]]$aic}
  if (crit=='icl'){ id_max = which.max(icl); crit_max = RES[[id_max]]$icl}
  #if (crit=='fisher'){ id_max = which.max(diff(fish)); crit_max = RES[[id_max]]$fish}
  res = RES[[id_max]]
  res$init_call = res0$call
  res$call = call
  res$crit = crit
  res$l1 = l1[id_max]
  res$l2 = l2
  if (disp){ 
    if (length(l1)>1) {
      cat('The best sparse model is with l1 penalty =',l1[id_max],'(',crit,'=',crit_max,')\n')
    } else { 
      cat('The sparse model has a l1 penalty =',l1[id_max],'(',crit,'=',crit_max,')\n')
    }
  }
  class(res)='bfem'
  res
}
