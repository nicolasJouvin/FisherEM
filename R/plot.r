

## Plot ##
#
# Author: 	C.Bouveyron & C.Brunet
# Institute:	SAMM - Paris 1 / Evry
# Date: 	February 2011
# License:	GPL v2
# 
# ***************************************************************
# plot function  : loglik + projection on the 2 first discriminative axes
# use: plot(obj) : - obj : object from fisherEM()  
# ***************************************************************
plot.fem <- function(x,loglik=F,...){	

Y = x$data
plot(as.data.frame(as.matrix(Y) %*% x$V[,1:2]),col=max.col(x$P),main='Estimated latent space with FisherEM',xlab='Discriminative axis 1',ylab='Discriminative axis 2',pch= max.col(x$P))
if (loglik){ x11()
	     plot(x$loglik,xlab='iterations',ylab='log-likelihood in the observed space',col=2,pch=20,cex=0.5,type="l")} 
	}

	