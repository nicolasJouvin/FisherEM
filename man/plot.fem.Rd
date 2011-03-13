\name{plot.fem}
\alias{plot.fem}

\title{
Plots the fem results
}
\description{
Plots the projection of clustered iris data into the estimated latent discriminative subspace with Fisher-EM 
or the evolution of associated log-likelihood
}
\usage{
plot.fem(x, loglik = F,...)
}

\arguments{
  \item{x}{
object from fem(), i.e, the clustering results obtained with the Fisher-EM algorithm
}
  \item{loglik}{
If loglik=F the function plots the projection of clustered iris data into the estimated latent discriminative subspace with Fisher-EM. 

If loglik=T the function plots the evolution of associated log-likelihood
}
\item{...}{
	Arguments to be passed to methods.
}
}




\author{
J.Loisel  T.Souvannarath  M.Tchebanenko  C.Bouveyron  C.Brunet
}

\examples{
Y<-iris[,1:4]
x<-fem(Y,3,model='AkB')
plot.fem(x,loglik=FALSE)
}