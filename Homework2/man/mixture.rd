\name{mixture}
\alias{mixture}
\title{Mixture Model Parameter Estimate} 
\description{
Estimate Parameters in two mixed normal distributions using Newton's method or EM algorithm.
}
\usage{
mixture(y, method=c("newton","EM"), maxit = NULL, tol = 1e-08, param0 = NULL)
}
\arguments{
\item{y}{Samples from the mixed normal distribution.}
\item{method}{Choose estimation method. Both newton or EM are valid.}
\item{maxit}{maximum number of iterations. If NULL, 100 max iterations for newton and 500 max iterations for EM.}
\item{tol}{Tolerance when the algorithm converges}
\item{param0}{Inital parameters. If NULL, the function will guess initial parameters according to y.}
}
\details{
Newton's method is not quite reliable. Both method are highly sensitive with the initial parameters. 
}
\value{
A list with elements mle containing the vector of maximum likelihood estimates and stderr containing the vector of corresponding asymptotic standard errors for the MLEs.
}

