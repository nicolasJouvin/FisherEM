bfem.init.prms <- function (Y,U,tau, model,method) {
  # Frequentist Mstep to initialize bfem.main and bem.sparse
  
  # Initialization
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  K = ncol(tau)
  U = as.matrix(U)
  d = ncol(U)
  mu   = matrix(NA,K,d)
  m   = matrix(NA,K,p)
  prop = rep(c(NA),1,K)
  D = array(0,c(K,d,d))
  b = rep(NA,K)
  
  # Projection
  X = as.matrix(Y) %*% as.matrix(U)
  
  # Estimation
  test = 0
  bk = 0
  for (k in 1:K){
    nk  = sum(tau[,k])
    if (nk ==0) stop("some classes become empty\n",call.=FALSE)
    # Prior Probability
    prop[k] = nk / n
    # Mean in the latent space
    mu[k,]  = colSums((tau[,k]*matrix(1,n,d))* X) / nk
    # Observed space
    # Needed for regular EM in DLM : \hat{m}_k = U %*% \mu_k
    m[k,]  = mu[k,] %*% t(U)
    YY  = as.matrix(Y - t(m[k,]*matrix(1,p,n)))
    if (nk < 1) denk = 1 else denk = (nk-1)
    YYk = tau[,k] * matrix(1,n,p) * YY
    # Estimation of Delta k amongst 12 submodels
    if (model %in% c('DkBk', 'DkB', 'DBk', 'DB')) {
      D[k,,] = crossprod(YYk%*%U, YY%*%U) / denk
    }
    else if (model %in% c('AkjBk', 'AkjB', 'AjBk', 'AjB')) {
      D[k,,] = diag(diag(crossprod(YYk %*%U, YY%*%U) / denk), d)
    }
    else if (model %in% c('AkBk', 'AkB', 'ABk', 'AB')) {
      D[k,,] = diag(rep(Trace(crossprod(YYk %*%U, YY%*%U) / denk)/d,d), d)
    }
    
    
    # compute b[k]
    if (model %in% c('DkB','DB','AkjB', 'AjB', 'AkB','AB')) {
      bk =  bk + (1 / n) * (sum(YYk*YY) - Trace(D[k,,]))
    } else {
      b[k] = (sum(YYk*YY) - Trace(D[k,,])) / (denk * (p-d))
    }
    
  }
  
  # Compute Sigma after loop for homoscedastic models
  Sigma = matrix(0, d, d) 
  if (model %in% c('DBk', 'DB', 'AjBk', 'AjB', 'AB', 'ABk')) {
    for (k in 1:K) {
      nk = sum(tau[,k])
      if (nk < 1) denk = 1 else denk = nk
      Sigma = Sigma + (denk / n) * D[k,,]
    }
    for (k in 1:K) {
      D[k,,] = Sigma
    }
  }
  
  # compute beta for '*B' models
  if (model %in% c('DkB','DB','AkjB','AkB','AjB','AB')) {
    bk = bk / (p-d)
    b = rep(bk, K)
  }
  
  # avoid numerical
  b[b<=0] = 1e-3
  for (k in 1:K) if (Trace(D[k,,]<1e-3)!=0) test = test+1
  
  prms = list(K=K,p=p,d=d,mean=mu,my=m,prop=prop,D=D,b=b,model=model,method=method,V=U,test=test)
}