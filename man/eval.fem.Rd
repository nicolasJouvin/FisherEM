\name{eval.fem}
\alias{eval.fem}

\title{
Computation of clustering accuracy when labels are known
}
\description{
Computes clustering accuracy when labels are known
}
\usage{
eval.fem(obj, cls, disp = 1)
}

\arguments{
  \item{obj}{
object from fem(), i.e, the clustering results obtained with the Fisher-EM algorithm
}
  \item{cls}{
vector of the class of each observation
}
  \item{disp}{
if disp is equal to integer value different than 0, the function displays us
"the Clustering accuracy rate" and "the Confusion matrix"
}
}

\value{
\item{lbl}{vector of the class of each observation}
\item{tx}{the Clustering accuracy rate}
}

\author{
J.Loisel  T.Souvannarath  M.Tchebanenko  C.Bouveyron  C.Brunet
}

\examples{
Y<-iris[,1:4]
obj<-fem(Y,3,model='AkB')
cls=as.numeric(iris[,5])
res=eval.fem(obj,cls,disp=1)
}


