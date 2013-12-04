\name{mixture}
\alias{mixture}

\title{
Estimate Unknown parameters from a mixture of 2 normal distribution
}
\description{
This package gives a function to estimate the 5 unknown parameters from a mixture of 2 normal distributions}
\usage{
mixture(y, method, maxit = NULL, tol = 1e-08, param0 = NULL)
}

\arguments{
  \item{y}{
dataset sampling from a mixture of 2 normal distributions}
  \item{method}{
"EM" stands for EM algorithm, "newton" stands for Newton's method'
}
  \item{maxit}{
maximum iteration numbers. By default, EM is 500, Newton is 100
}
  \item{tol}{
level of tolerance. The default level is 1e-8
}
  \item{param0}{
initialization parameter
}
}
\details{
You need to set the desired method when using this function. Newton's method requires a good initialization parameter. Therefore, it would be better to use default set
}
\value{
  \item{mle}{max likelihood parameters}      
  \item{stderr}{asymptotic stderr}
  
  }
\references{
PH140.778, Professor Roger Peng
}
\author{
Lu Li}