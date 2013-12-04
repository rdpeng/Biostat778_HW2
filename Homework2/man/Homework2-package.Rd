\name{Homework2-package}
\alias{Homework2-package}
\alias{Homework2}
\docType{package}
\title{
Mixture of two normal distributions
}
\description{
Approximate mixed normal distribution using Newton or EM algorithm
}
\details{
\tabular{ll}{
Package: \tab Homework2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-04\cr
License: \tab GPL \cr
}
}
\author{
Yuan He

Maintainer: Yuan He <yhe23@jhu.eud>

}
\references{
https://www.cs.duke.edu/courses/spring04/cps196.1/handouts/EM/tomasiEM.pdf
}
\keyword{ EM }
\examples{
y = c(rnorm(100,10,10),rnorm(19,50,40))
mixture(y,method="EM")
}
