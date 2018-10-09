\name{fem.ari}
\alias{fem.ari}
\title{
Adjusted Rand index
}
\description{
The function computes the adjusted Rand index (ARI) which allows to compare two clustering partitions.}
\usage{
fem.ari(x,y)
}

\arguments{
  \item{x}{
A 'fem' object containing the first partition to compare.
}
  \item{y}{
The second partition to compare (as vector).
}
}

\value{
  \item{ari}{The value of the ARI.}
}

\seealso{fem, sfem, plot.fem, summary.fem}
\examples{
data(iris)
res = fem(iris[,-5],K=3,model='DkBk',method='reg')
res
plot(res)
fem.ari(res,as.numeric(iris[,5]))
}

