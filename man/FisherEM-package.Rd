\name{FisherEM-package}
\alias{FisherEM-package}
\alias{FisherEM}
\docType{package}
\title{
Model-based Clustering in the Fisher Discriminative Subspace
}
\description{
Model-Based Clustering in the Fisher Discriminative Subspace provides 
a low-dimensional discriminative representation of the clustered data.
To find a parsimonious and discriminative fit for the data this method 
uses discriminative latent model (DLM).The Fisher-EM algorithm estimates 
the parameters of DLM models in order to cluster and visualize the clustered data.
}
\details{
\tabular{ll}{
Package: \tab FisherEM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-03-13\cr
License: \tab GPL-2 \cr
LazyLoad: \tab yes\cr
}
This package is used to make simultaneous model-based clustering with high 
dimensional data using discriminative latent mixture (DLM) model with 
constraining parameters within and between groups. Twelve parsimonious 
DLM models were set up to provide various situations. The fem() function
computes a discriminant subspace alternatively with Gaussian Mixture Model 
and using the Fisher-EM algorithm which is divided into three steps due to the 
nature of the model. The F step estimates the orientation matrix 
of the discriminative latent space conditionally to the posterior 
probalities. The E step computes the posterior probabilities of 
each observations. The M step updates the parameters of the 
Mixture model (mean, prior probabilities).
}
\author{
J.Loisel  T.Souvannarath  M.Tchebanenko  C.Bouveyron  C.Brunet

}
\references{
C.Bouveyron and C.Brunet, Simultaneous model-based clustering and visualization in the Fisher discriminative subspace, Statistics and Computing, to appear, 2011

}




