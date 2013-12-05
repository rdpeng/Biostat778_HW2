\name{mixture}
\alias{mixture}
\title{Mixture models function}
\usage{
  mixture(y, method = "EM", maxit = 10000, tol = 1e-08,
    param0 = NULL)
}
\arguments{
  \item{y}{data vector}

  \item{method}{either 'EM' or 'newton'}

  \item{maxit}{number of iterations for either algorithm}

  \item{tol}{how much any parameter can differ from the
  previous fitted parameter. If any parameter differs by
  more than this much, we keep searching}

  \item{param0}{A starting value for the MLE. If none are
  provided, we calculate a nieve estimate using the
  quantiles of y.}
}
\description{
  \code{mixture} Fits a mixture model with either EM or
  newton
}

