\name{plot.fem}
\alias{plot.fem}

\title{
The plot.fem function
}
\description{
This function enables to plot the clustered data in the two axes the most discriminative fitted by the Fisher-EM algorithm.
}
\usage{
\method{plot}{fem}(x, Y, loglik=F, ...)
}

\arguments{
  \item{x}{the fem object
}
  \item{Y}{
the data matrix containing the individuals in lines and the variables in columns
}
  \item{loglik}{
enables to plot the loglikelihood values according to the number of iterations if loglik=T. By default, loglik=F.
}
  \item{...}{
some additional options for the plot.
}
}


\seealso{fem, evalEM}

\examples{
data(iris)
res = fem(iris[,-5],3,model='AkB')
plot.fem(res,iris[,-5])
}