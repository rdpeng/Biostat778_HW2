\name{mixture}
\alias{mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of Gaussian mixture model parameters}
\description{
Estimates the mixing parameter, means, and variances for the mixture of 2 normal densities:
\eqn{f(y_i) = \lambda \phi(y_i; \mu_1,\sigma_1^2) + (1-\lambda)\phi(y_i;\mu_2,\sigma_2^2)}
}
\usage{
mixture(y, method, maxit = NULL, tol = 1e-08, param0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  Vector of observations from the mixture model
}
  \item{method}{
  Choice of method used to fit the mixture model. Choices include “newton” for Newton's method or “EM” for the Expectation-Maximization algorithm
}
  \item{maxit}{
  Maximum number of iterations to be performed by the chosen method. Default is NULL which indicates a maximum of 100 iterations for Newton's method and 500 for the EM algorithm.
}
  \item{tol}{
  Controls the tolerance for convergence (default is 1e-8)
}
  \item{param0}{
  Starting value for Newton's method or the EM algorithm. Default is  NULL, in which case the starting value is automatically determined by the function.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{coefficients}{A vector of the estimated beta coefficients}
\item{vcov}{The variance-covariance matrix of the estimated beta coefficients}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Leslie Myint
}
\examples{
set.seed(1)
y1 = rnorm(100,5,1)
y2 = rnorm(200,9,2)
y = c(y1,y2)
fitnewt = mixture(y,"newton")
fitem = mixture(y,"EM")
}