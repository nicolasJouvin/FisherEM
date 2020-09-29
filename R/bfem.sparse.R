bfem.sparse <- function(Y,K,maxit.sparse,eps.sparse,maxit.ve, eps.ve, Tinit, 
                        lambda, nu, emp.bayes,
                        model,method='reg',l1,nbit,l2){
  colNames = colnames(Y)
  if (length(l1)!=1 | l1>1) 
    stop('The l1 penalty term is a single figure comprises between 0 and 1')
  
  ve.tol = eps.ve
  ve.max.iter = maxit.ve
  
  # Initialization
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  d = min((K-1),(p-1))
  
  # New objects
  Lobs = rep(c(-Inf),1,(maxit.sparse+1))
  	
  # Initialization of tau
  tau       = Tinit
  logtau    = log(tau)
  
  # Init of sparse U
  U         = fstep.sparse(Y,tau,l1,nbit,l2)
  
  # Init Sigmak and Beta with frequentist M-step
  prms.fem  =  bfem.init.prms(Y,U,tau,model=model,method=method)
  prms      = transform_param_fem_to_bfem(prms.fem)
  
  # Inital VE-step: 1 iteration of fixed point for var parameter of q(\mu_k)
  Varcovk   = updateVarcovk(prms, lambda, tau)
  Varmeank  = updateVarmeank(Y, U, prms, nu, tau, Varcovk)
  
  elbo      = bfem.elbo(Y, U, prms, nu, lambda, logtau, Varmeank, Varcovk)
  
  # ================ Main loop BFEM ===============================
  Lobs = rep(-Inf , 1, maxit.sparse + 1)
  Lobs[1] = elbo
  Linf_new  = Lobs[1]
  for (i in 1:maxit.sparse){
    # The three main steps F, VE and M
    U         = fstep.sparse(Y,tau,l1,nbit,l2)
    
    # Variational E-step (don't change tau yet, only q(\mu_k))
    ve_step   = bfem.vestep(Y, U, prms, nu, lambda, logtau, 
                          Varmeank=NULL, Varcovk=NULL, ve.tol, ve.max.iter)
    Varmeank = ve_step$Varmeank
    Varcovk = ve_step$Varcovk   
    
    # Mstep
    prms = bfem.mstep(Y, U, tau, Varmeank, Varcovk, model)
    if(is.null(prms$Sk)) return(list(k=K, Tinit=Tinit, loglik=-Inf, 
                                     icl = -Inf, bic = -Inf, aic = -Inf, 
                                     model = model, K=K, cls = max.col(tau), P = tau,
                                     error=1))
    
    # Update of tau and logtau after M-step
    logtau = updatelogTau(Y, U, prms, Varmeank, Varcovk) 
    tau = exp(logtau)
    
    # Hyper-params 
    if (emp.bayes) {
      nu = updateNu(Varmeank)
      lambda = updateLambda(Varmeank, Varcovk, nu)
    }
    
    # Stop criterion
    elbo = bfem.elbo(Y, U, prms, nu, lambda, logtau, Varmeank, Varcovk)
    
    if (is.na(elbo)) {
      warning('Likelihood becomes NA (probably emptied cluster).')
      break
    }
    Lobs[i+1] = elbo

    if (i>=2){
      acc = (Lobs[i+1] - Lobs[i]) / (Lobs[i] - Lobs[i-1])
      Linf_old = Linf_new
      Linf_new <- try( Lobs[i] + 1/(1 - acc) * (Lobs[i+1] - Lobs[i]))
      if (is.na(Linf_new)){warnings('acc=', acc);break}
      if (abs(Linf_new - Linf_old) < eps.sparse) {break}
    }
  }
  Lobs = Lobs[Lobs != -Inf]
  
  # Returning the results
  cls = apply(tau, 1, which.max) 
  param = list(PI = prms$PI, Sigmak = prms$Sk$Sigmak, Beta = prms$Sk$Beta)
  var_param = list(logtau = logtau,  Varmeank = Varmeank, Varcovk = Varcovk)
  hypers = list(lambda = lambda, nu = nu)
  proj = Y %*% U
  rownames(U) = colNames
  final_elbo = utils::tail(Lobs, 1)
  crit = bfem.criteria(Y, final_elbo, U, prms, tau, model, hypers, var_param, 
                       method = 'sparse')
  return(list(K = K, Tinit = Tinit, d = d, cls = cls, P = tau,
              elbos = Lobs, loglik=final_elbo, n_ite = i,
              U = U, param = param, var_param = var_param, proj = proj,
              model = model, aic=crit$aic, bic=crit$bic, icl=crit$icl,
              comp=crit$comp,
              hypers = hypers,
              call = call)
         )
}

