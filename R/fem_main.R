library('MASS')
fem_main <-
function(Y,K,init,maxit,eps,Tinit,model,kernel,method,graph,Hess){

  # Initialization
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  d = min((K-1),(p-1))
  # 	lambda = c() # not allowed in the global function
  # 
  # New objects
  Lobs = rep(c(-Inf),1,(maxit+1))
  # 	
  # Initialization of T
  if (init=='user'){ T = Tinit}
  if (init=='kmeans'){
    T   = matrix(0,n,K)
    ind = kmeans(Y,K)$cluster
    for (i in 1:n){ T[i,ind[i]] = 1 }
  }
  if (init=='random'){
#     set.seed(100)
    T   = t(rmultinom(n,1,c(rep(1/K,K))))
  }
#   if (init=='miniEM'){
#     source('EM.r')
#     TT = array(NA,c(10,n,K))
#     Lik = rep(NA,10)
#     for (nb in 1:10){
#       res.em   = EM(Y,K,model='full',maxit=10)
#       Lik[nb]  = res.em$lik
#       TT[nb,,] = res.em$P
#     }
#     T   = TT[which.max(Lik),,]
#   }
  V = switch(method,
             'SVD'= fstep.fisher(Y,T,kernel),
             'GS'= fstep.GramSc(Y,T,kernel),
             'REG'  = fstep.qiao(Y,T,kernel))
  prms      = mstep(Y,V,T,model=model,method=method)
  res.estep = estep(prms,Y,V)
  T         = res.estep$T
  Lobs[1]   = res.estep$loglik
  
  # Main loop
  Linf_new  = Lobs[1]
  for (i in 1:maxit){
#     cat('*')
    # sum(diag())
    # The three main steps F, M and E
    V = switch(method,
               'SVD'= fstep.fisher(Y,T,kernel),
               'GS'= fstep.GramSc(Y,T,kernel),
               'REG'  = fstep.qiao(Y,T,kernel))
    prms      = mstep(Y,V,T,model=model,method=method)
    res.estep = estep(prms,Y,V)
    T         = res.estep$T
    Lobs[i+1] = res.estep$loglik
    
    # Stop criterion
    if (i>=2){
      acc      = (Lobs[i+1] - Lobs[i]) / (Lobs[i] - Lobs[i-1])
      Linf_old = Linf_new
      Linf_new <- try( Lobs[i] + 1/(1-acc) * (Lobs[i+1] - Lobs[i]))
      if (abs(Linf_new - Linf_old) < eps) {break}
    }
    
  }
#   cat('\n')
  
  # Compute Hessian to verify the convergence
  H = NA
  if (Hess){
    theta = c(prms$prop[1:(K-1)],prms$D[,1,1],prms$D[,p,p],as.vector(t(prms$my)))
    Hess  = fdHess(theta,fem.loglik,K,Y,V)
    H     = Hess$Hessian
    H     = diag(1,ncol(H)) - H
    H     = eigen(H)$val }
  
  # Graphic option
  if (graph){
    par(mfrow=c(1,2))
    plot(as.data.frame(as.matrix(Y) %*% V[,1:2]),col=max.col(T),xlab='axis 1',ylab='axis 2',pch=20)
    plot(Lobs[1:i],xlab='iterations',ylab='Likelihood',col=2,pch=20)
  }
  
  # Returning the results
  cls  = max.col(T)
  crit = criteria(Lobs[(i+1)],T,prms,n); # BIC
  res  = list(cls=cls,P=T,prms=prms,U=V,aic=crit$aic,bic=crit$bic,icl=crit$icl,loglik=Lobs[2:(i+1)],ll=Lobs[i+1],Hess=H,method=method)
  
  res
}

