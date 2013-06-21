\name{summary.fem}
\alias{summary.fem}

\title{
The summary function for 'fem' objects.
}
\description{
This function summarizes 'fem' objects. It in particular indicates which DLM model has been chosen and displays the loading matrix 'U' if the original dimension is smaller than 10.
}
\usage{
 \method{summary}{fem}(object,...)
}

\arguments{
  \item{object}{
The fem object.
}
  \item{...}{
Additional options to pass to the summary function.
}
}

\seealso{fem, sfem, fem.ari, plot.fem}

\examples{
data(iris)
res = fem(iris[,-5],K=2:5,model='AkB')
summary(res)
plot(res)
fem.ari(res,as.numeric(iris[,5]))
}
