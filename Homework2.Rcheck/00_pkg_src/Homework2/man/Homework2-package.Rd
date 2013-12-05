\name{Homework2-package}
\alias{Homework2-package}
\alias{Homework2}
\docType{package}
\title{
Homework 2 for Advanced Computing
}
\description{
The package contains a function, 'mixture', which applies Newton and EM algorithm to estimate Gaussian mixture model parameters.
}
\details{
\tabular{ll}{
Package: \tab Homework2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-04\cr
License: \tab GPL-2\cr
}
mixture(y, method, maxit = NULL, tol = 1e-08, param0 = NULL)
input data y, and choose method.
}
\author{
Detian Deng

Maintainer: Detian Deng <ddeng@jhsph.edu>

}
\references{
class notes 
}

\examples{
w = 0.6 
m1 = 5; m2 = -5; v1 = 2; v2 = 3

W = rbinom(1000,1,prob=w)
N1 = rnorm(1000,m1,sqrt(v1)) 
N2 = rnorm(1000,m2,sqrt(v2))
Y = W*N1+(1-W)*N2
rm(w,m1,m2,v1,v2,W,N1,N2)

mixture(y=Y,method="Newton",maxit = 100)
mixture(y=Y,method="EM",maxit = 300)
}
