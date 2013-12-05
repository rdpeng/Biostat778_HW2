\name{mixture}
\alias{mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate mle for gaussian mixture
%%  ~~function to do ... ~~
}
\description{This function generates the mle for gaussian mixture with a set of ovserved data.There are two options.One is Em algorithm and the other is Newton method.Also, for each method,standard error is given.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mixture(y, method, maxit = NULL, tol = 1e-08, param0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The ovserved data
%%     ~~Describe \code{y} here~~
}
  \item{method}{The method to use
%%     ~~Describe \code{method} here~~
}
  \item{maxit}{the maximum number of iterations
%%     ~~Describe \code{maxit} here~~
}
  \item{tol}{the tolerance of the iterations
%%     ~~Describe \code{tol} here~~
}
  \item{param0}{the starting point of the iterations
%%     ~~Describe \code{param0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list of 2. One is mle and the other is standard error.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Jerryppa
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	a1<-rbinom(1000,1,0.5)
a2<-vector(length=1000)
for (i in 1:1000){
  if(a1[i]==0){
    a2[i]=rnorm(1,mean=1,sd=1)
    
  }
  else{a2[i]=rnorm(1,mean=2,sd=2)}
  
}

mixture(a2,method="EM",maxit=100)
	
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, method, maxit = NULL, tol = 1e-08, param0 = NULL) 
{
    mle = vector(length = 5)
    stderr = vector(length = 5)
    if (is.null(param0) == TRUE) {
        lambda1 = 0.5
        lambda2 = 0.5
        mu1 = 10.5
        sigma1 = 59
        mu2 = 20.5
        sigma2 = 249
    }
    else {
        lambda1 = param0[1]
        mu1 = param0[2]
        mu2 = param0[3]
        sigma1 = param0[4]
        sigma2 = param0[5]
    }
    m <- length(y)
    if (method == "EM") {
        it = 0
        if (is.null(maxit) == TRUE) {
            maxit = 500
        }
        T <- matrix(rep(0, m * 2), nrow = m, ncol = 2)
        for (i in 1:maxit) {
            ll = sum(log(lambda1 * dnorm(y, mean = mu1, sd = sqrt(sigma1)) + 
                lambda2 * dnorm(y, mean = mu2, sd = sqrt(sigma2))))
            f1 = dnorm(y, mean = mu1, sd = sqrt(sigma1))
            f2 = dnorm(y, mean = mu2, sd = sqrt(sigma2))
            T[, 1] = lambda1 * f1/((lambda1 * f1) + (lambda2 * 
                f2))
            T[, 2] = lambda2 * f2/((lambda1 * f1) + (lambda2 * 
                f2))
            lambda1 = mean(T[, 1])
            lambda2 = mean(T[, 2])
            mu1 = sum(T[1:m, 1] * y[1:m])/sum(T[1:m, 1])
            mu2 = sum(T[1:m, 2] * y[1:m])/sum(T[1:m, 2])
            sigma1 = sum(T[1:m, 1] * (y[1:m] - mu1)^2)/sum(T[1:m, 
                1])
            sigma2 = sum(T[1:m, 2] * (y[1:m] - mu2)^2)/sum(T[1:m, 
                2])
            lln = sum(log(lambda1 * dnorm(y, mean = mu1, sd = sqrt(sigma1)) + 
                lambda2 * dnorm(y, mean = mu2, sd = sqrt(sigma2))))
            if (abs(lln - ll) <= tol) {
                break
            }
            it = it + 1
        }
        dl = T[, 1]/lambda1 - T[, 2]/lambda2
        dmu1 = T[, 1] * (y - mu1)/sigma1
        dmu2 = T[, 2] * (y - mu2)/sigma2
        ds1 = 0.5 * T[, 1] * ((y - mu1)^2 - sigma1)/sigma1^2
        ds2 = 0.5 * T[, 2] * ((y - mu2)^2 - sigma2)/sigma2^2
        temp = rbind(dl, dmu1, dmu2, ds1, ds2)
        ss = matrix(rep(0, 25), nrow = 5, ncol = 5)
        for (i in 1:m) {
            ss = ss + temp[, i] \%*\% t(temp[, i])
        }
        va = sqrt(diag(solve(ss)))
        mle = c(lambda1, mu1, mu2, sigma1, sigma2)
        stderr = va
    }
    else if (method == "newton") {
        if (is.null(maxit) == TRUE) {
            maxit = 100
        }
        x = c(lambda1, mu1, mu2, sigma1, sigma2)
        it = 0
        for (i in 1:maxit) {
            ll = sum(log(lambda1 * dnorm(y, mean = mu1, sd = sqrt(sigma1)) + 
                lambda2 * dnorm(y, mean = mu2, sd = sqrt(sigma2))))
            gd = der(x, y, m)$gd
            hes = der(x, y, m)$hes
            x = x - solve(hes) \%*\% gd
            lambda1 = x[1]
            mu1 = x[2]
            mu2 = x[3]
            sigma1 = x[4]
            sigma2 = x[5]
            lln = sum(log(lambda1 * dnorm(y, mean = mu1, sd = sqrt(sigma1)) + 
                lambda2 * dnorm(y, mean = mu2, sd = sqrt(sigma2))))
            if (abs(lln - ll) <= tol) {
                break
            }
            it = it + 1
        }
        lambda1 = x[1]
        mu1 = x[2]
        mu2 = x[3]
        sigma1 = x[4]
        sigma2 = x[5]
        va = der(x, y, m)$vr
        mle = c(lambda1, mu1, mu2, sigma1, sigma2)
        stderr = va
    }
    mle = as.matrix(mle)
    stderr = as.matrix(stderr)
    rownames(mle) = c("lambda", "mu1", "mu2", "sigma1", "sigma2")
    rownames(stderr) = c("lambda", "mu1", "mu2", "sigma1", "sigma2")
    list(mle = mle, stderr = stderr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
