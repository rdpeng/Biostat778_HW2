\name{Homework2-package}
\alias{Homework2-package}
\alias{Homework2}
\docType{package}
\title{Generate mle for gaussian mixture
What the package does (short line)
~~ package title ~~
}
\description{This function generates the mle for gaussian mixture with a set of ovserved data.There are two options.One is Em algorithm and the other is Newton method.Also, for each method,standard error is given.
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab Homework2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-04\cr
License: \tab What license is it under?\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{Jerryppa
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
~~ directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{a1<-rbinom(1000,1,0.5)
a2<-vector(length=1000)
for (i in 1:1000){
  if(a1[i]==0){
    a2[i]=rnorm(1,mean=1,sd=1)
    
  }
  else{a2[i]=rnorm(1,mean=2,sd=2)}
  
}

mixture(a2,method="EM",maxit=100)
}
