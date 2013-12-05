\name{mixture}
\alias{mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimating Gaussian mixture model parameters
%%  ~~function to do ... ~~
}
\description{
Apply Newton and EM algorithm to estimate Gaussian mixture model parameters.
}
\usage{
mixture(y, method, maxit = NULL, tol = 1e-08, param0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Data in a vector or a data frame.
}
  \item{method}{
"Newton" or "EM".
}
  \item{maxit}{
Maximum number of iteration allowed.
}
  \item{tol}{
Convergence tolerance. 
}
  \item{param0}{
Starting values of parameters. Default to NULL.
}
}

\value{
\item{mle}{Maximul Likelihood Estimate of parameters.}
\item{stderr}{Standard error of the MLEs}
%% ...
}

\author{
Detian Deng
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

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

