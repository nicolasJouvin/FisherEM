\name{evalEM}
\alias{evalEM}
\title{
The evalEM function
}
\description{
evalEM provides a contingency table and a correct classification rate from the clustering results of the Fisher-EM algorithm and true labels.
}
\usage{
evalEM(cls, obj, disp = 1)
}

\arguments{
  \item{cls}{
is a numerical vector containing the true labels for each individual
}
  \item{obj}{
the object obtained by the fem() function
}
  \item{disp}{
if disp=1 then the contingency table is displayed.
}
}

\value{
  \item{lbl}{the permuted groups}
  \item{tx}{the correct classification rate}
}

\seealso{fem, plot.fem}
\examples{
## The function is currently defined as

data(iris)
cls = as.numeric(iris[,5]); # the true numerical labels
res = fem(iris[,-5],3,model='AkB');
evalEM(cls,res)
}

