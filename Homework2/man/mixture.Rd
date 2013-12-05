\name{mixture}
\alias{mixture}
\title{
Mixture
}
\description{
Calculate MLEs using NR and EM
}
\usage{
mixture(y, method, maxit = NULL, tol = 1e-08, param0 = NULL)
}

\arguments{
  \item{y}{inputing data}
  \item{method}{choosing either "Newton" or "EM"}
  \item{maxit}{maximum number of iteration}
  \item{tol}{tolerence}
  \item{param0}{starting value}
}