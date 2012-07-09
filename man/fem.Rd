\name{fem}
\alias{fem}

\title{
The Fisher-EM algorithm
}
\description{
Fem is a subspace clustering method. It is based on the
Gaussian Mixture Model and on the idea that the data lives in
a common and low dimensional subspace. An EM-like algorithm estimates both 
the discriminative subspace and the parameters of the mixture model.
}
\usage{
fem(Y,K,init='random',maxit=100,eps=1e-6,Tinit=c(),model='AkjBk',
      kernel='',graph=F,Hess=F,method='SVD',l1=0.3,l2=0,nbit=2)
}
\arguments{
  \item{Y}{
contains the data matrix (without NAs)
}
  \item{K}{
the number of clusters
}
  \item{init}{
the kind of initialization of the Fisher-EM algorithm. There are 3 options: "random" for a randomized initialization, "kmeans" for an initialization by the traditional kmeans algorithm
or "user" for a chosen initialization for which the parameter Tinit needs to be completed.
}
  \item{maxit}{
the maximum number of iterations before the stop of the Fisher-EM algorithm. 
}
  \item{eps}{
the threshold for the stop of the Fisher-EM algorithm.
}
  \item{Tinit}{
a n x K matrix which contains posterior probabilities: each line corresponds to an individual 
}
  \item{model}{ the kind of Discriminative Latent Mixture model. There are 12 different models: "DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk", "AjBk",
"AjB", "ABk", "AB".  The option "all" executes the Fisher-EM algorithm on the 12 DLM models and select the best model according to the maximum value obtained by 
BIC criterion. 
}
  \item{kernel}{ enables to deal with the n < p problem. By default, no kernel ' ' is used. But the user has also the choice between 3 options for the kernel: 'linear', 
'sigmoid' or 'rbf'.
}
  \item{graph}{
if TRUE, the clustered data are plotted on the 2 first discriminative axes fitted by the Fisher-EM algorithm.
}
  \item{Hess}{
if TRUE, the Hessian matrix is computed.
}
  \item{method}{
corresponds to 3 different methods ('SVD', 'REG', 'GS') for the fitting of the projection matrix associated to the discriminative subspace. The 'Fisher' method 
is used by default. If the option method='sparse' is used, then the loadings of the projection matrix are sparse. The parameter l1 enables to change the level of sparsity.
}
  \item{l1}{
is a l1 penalty term: the figure is comprised between 0.1 (very sparse loadings of the projection matrix) and 1 (no sparsity). This option has to be used with 
method='sparse'.
}
  \item{l2}{
is a l2 penalty term.
}
  \item{nbit}{
number of iterations for fitting the sparse loadings at each update of the projection matrix.
}
}


\value{
A list is returned: 
  \item{cls}{the group membership of each individual estimated by the Fisher-EM algorithm}
  \item{P}{the posterior probabilities of each individual for each group}
  \item{prms}{lists of parameters of the mixture model fitted by the Fisher-EM algorithm.}
  \item{U}{projection matrix}
  \item{aic}{Akaike criterion}
  \item{bic}{Bayesian Information criterion}
  \item{icl}{Integrated Completed Likelihood criterion}
  \item{loglik}{log-likelihood values computed at each iteration of the FEM algorithm}
  \item{ll}{the log-likelihood value obtained at the last iteration of the FEM algorithm}
  \item{Hess}{the Hessian matrix if Hess = TRUE}
  \item{method}{Method used}
}
\references{Charles Bouveyron and Camille Brunet, Simultaneous model-based clustering and visualization in the Fisher discriminative subspace,
        Statistics and Computing, 22(1), 301-324.

	    Charles Bouveyron, Camille Brunet (2012), "Discriminative variable selection for clustering with the sparse Fisher-EM algorithm",
        preprint Hal n-00685183.

}
\author{
Charles Bouveyron and Camille Brunet
}

\seealso{plot.fem, evalEM}

\examples{
data(iris)
res1 = fem(iris[,-5],3,model='AkB')
res1$U # print the loadings of the projection matrix 

## For a sparse case:
res2 = fem(iris[,-5],3,model='AkB',method='sparse', l1=0.2)
res2$U #print the loadings of the projection matrix
}

%%\keyword{ gmm }
%%\keyword{ subspace clustering }
%%\keyword{ high-dimensional data }
%%\keyword{ discriminative subspace }
%%\keyword{ variable selection }
