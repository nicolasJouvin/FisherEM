plot.fem <-
function(x,Y,loglik=F,...){
Y = as.matrix(Y)
p = ncol(Y)
n = nrow(Y)
	if (is.matrix(x$U)){
	    plot(as.data.frame(as.matrix(Y) %*% x$U[,1:2]),col=max.col(x$P),xlab='axis 1',ylab='axis 2',pch=20,...)} 
	else{
	      cat('\n', 'Since K=2, the data and the discriminative subspace have been projected on the 2 first PCs','\n')
		      if (ncol(Y)<nrow(Y))  Z = -eigen(cov(Y),symmetric=T)$vectors[,1:2]  
		      else {
			    z = eigen(Y%*%t(Y),symmetric=T)
			    Z = matrix(NA,p,2)
			    for (i in 1:2) Z[,i] = t(Y)%*%z$vectors[,i]/sqrt(n*z$values[i]) 
		      }

			  MU	= colMeans(Y)
			  proj= matrix(NA,2,p)
			  axU = matrix(x$U,p,1)%*%matrix(x$U,1,p)%*%matrix(10, p, 1)
			  proj= axU+matrix(MU,p,1)
			  Yproj = Y%*%Z
			  u 	= matrix(proj,1,p)%*%Z # projection sur les 2 pc
			  ybar 	= matrix(MU,1,p) %*% Z # proj des moyennes des cls sur les 2 pc
			  plot(Yproj,col=x$cls+1,xlab='comp.1',ylab='comp.2',pch=20)
			  pente=(u[i,2]-ybar[i,2])/(u[i,1]-ybar[i,1])
			  oo=u[i,2]-pente*u[i,1]
			  xb=(2*ybar[i,1]-sqrt(50^2/(pente^2+1)))/2
			  xa=(2*ybar[i,1]+sqrt(50^2/(pente^2+1)))/2
			  lines(c(xa,xb),oo+pente*c(xa,xb),col=1,type='l',lwd=2)
	}
	
	if (loglik){ x11()
		plot(x$loglik,xlab='iterations',ylab='log-likelihood in the observed space',col=2,pch=20,cex=0.5,...)
	} 
}

