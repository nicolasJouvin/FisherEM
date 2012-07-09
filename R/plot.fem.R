plot.fem <-
function(x,Y,loglik=F,...){	
	plot(as.data.frame(as.matrix(Y) %*% x$U[,1:2]),col=max.col(x$P),xlab='axis 1',ylab='axis 2',pch=20,...)
	if (loglik){ x11()
		plot(x$loglik,xlab='iterations',ylab='log-likelihood in the observed space',col=2,pch=20,cex=0.5,...)
	} 
}

