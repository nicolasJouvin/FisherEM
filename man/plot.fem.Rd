\name{plot.fem}
\alias{plot.fem}

\title{
The plot function for 'fem' objects.
}
\description{
This function plots different information about 'fem' objects such as model selection, log-likelihood evolution and visualization of the clustered data into the discriminative subspace fitted by the Fisher-EM algorithm.}
\usage{
 \method{plot}{fem}(x, frame=0, crit=c(),...)
}

\arguments{
  \item{x}{The fem object.
}
  \item{frame}{0: all plots; 1: selection of the number of groups; 2: log-likelihood; projection of the data into the discriminative subspace.
}
\item{crit}{
The model selection criterion to display. Default is the criterion used in the 'fem' function ('icl' by default).
}
  \item{...}{
Additional options to pass to the plot function.
}
}

\seealso{fem, sfem, fem.ari, summary.fem}

\examples{
data(iris)
res = fem(iris[,-5],K=2:5,model='AkB')
summary(res)
plot(res)
fem.ari(res,as.numeric(iris[,5]))
}
