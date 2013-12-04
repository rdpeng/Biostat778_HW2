\name{mixture}
\alias{mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MLE and stderr of mixture model parameters
}
\description{
calculates MLE and standard error for mixture model parameters using either EM or newton method.
}
\usage{
mixture(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
data from mixture distribution to be estimated;
}
\item{method}{
indicates which method to use, values can be "EM" or "newton";
}
\item{maxit}{
maximum number of iterations, default is NULL;
}
\item{tol}{
convergence criteria, default is 1e-08;
}
\item{param0}{
starting value, default is NULL.
}
}

\value{
a list is returned with MLE and standard error for lambda, mu1, mu2, sigma1, and sigma2.
}

\author{
Yuxin(Daisy) Zhu
}

\keyword{ mixture normal}
\keyword{  EM}% __ONLY ONE__ keyword per line
\keyword{newton}