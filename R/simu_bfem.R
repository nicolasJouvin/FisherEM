simu_bfem <- function(n, which = "section4.2", ...) {
  simu = switch(which,
                "Chang1983" = simu.Chang1983(n),
                "section4.2" = simu.section4.2(n, ...),
                "section4.3" = simu.section4.3(n, ...))
  return(simu)
}

# Chang 1983 dataset
simu.Chang1983 <- function(n = 300) {
  p = 15
  
  cte = matrix(rep(0.95 - 0.05* (1:p), n), nrow = n, byrow = T)
  Z = rbinom(n=n, size=1, prob = 0.2)
  
  mu = rep(0, p)
  Sigma = diag(1, p)
  f = c(rep(-0.9, 8), rep(0.5, 7))
  for(i in 1:p) {
    for (j in 1:p) {
      if (j!= i) Sigma[i,j] = -0.13 * f[i] * f[j] 
    }
  }
  X = MASS::mvrnorm(n=n, mu = mu, Sigma = Sigma)
  
  Y = 0.5 * cte + Z * cte + X
  return(list(Y=Y, cls = factor(Z)))
}

# High-dimensional setting of Section 4.2
simu.section4.2 <- function(n, p = 50, noise = 1) {
  d = 2
  K = 3
  S = array(0,c(d,d,K))
  S[,,2] = S[,,3] = S[,,1] = sigma * matrix(c(1.5,0.75, 0.75, 0.45))
  
  cst = 2.5
  mu = matrix(0, d, K)
  mu[,1] = c(0, 0)
  mu[,2] = cst * c(1, 0)
  mu[,3] = cst * c(2, 0)

  n1 = round(4*n/10)
  n2 = round(3*n/10)
  n3 = round(3*n/10)
  
  ind = sample(1:n)
  X = rbind(mvrnorm(n1,mu[,1],S[,,1]),mvrnorm(n2,mu[,2],S[,,2]),mvrnorm(n3,mu[,3],S[,,3]))
  X = X[ind,]
  A = mvrnorm(p, mu = rep(0,p),Sigma = diag(25, p, p))
  W = qr.Q(qr(A))
  B = mvrnorm(n,rep(0, p-d), diag(noise, p-d, p-d))
  y = cbind(X,B)
  Y = y %*% t(W)
  
  cls = c(rep(1,n1),rep(2,n2),rep(3,n3))[ind]
  
  return(list(Y = Y, W = W, X = X, cls = cls))  
}

# Signal-to-noise ratio of section 4.3
simu.section4.3 <- function(n, p = 150, snr = 3) {
  # inverse formula: snr = 10*log10(1.95 / beta)
  beta = 1.95 / (exp(snr * log(10)/10))
  simu.section4.2(n, p = p, noise = beta)
}