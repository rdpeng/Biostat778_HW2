\name{der}
\alias{der}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{To generate the gradient and Hessian of the likelihood function
%%  ~~function to do ... ~~
}
\description{To generate the gradient and Hessian of the likelihood function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
der(xt, y, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xt}{
%%     ~~Describe \code{xt} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{it contains several elements.One important element is gradient vector and another one is hessian matrix
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
	p0=c(0.6,10,20,60,250)
	y=c(10,15,20,20,15)
	l=length(y)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (xt, y, m) 
{
    gd = matrix(rep(0, 5), nrow = 1, ncol = 5)
    hes = matrix(rep(0, 25), nrow = 5, ncol = 5)
    nor = expression(log(lambda1/sqrt(2 * pi)/sqrt(sigma1) * 
        exp((-1) * (y - mu1)^2/(2 * sigma1)) + (1 - lambda1)/sqrt(2 * 
        pi)/sqrt(sigma2) * exp((-1) * (y - mu2)^2/(2 * sigma2))))
    gra = deriv3(nor, c("lambda1", "mu1", "mu2", "sigma1", "sigma2"))
    lambda1 = xt[1]
    mu1 = xt[2]
    mu2 = xt[3]
    sigma1 = xt[4]
    sigma2 = xt[5]
    Gra = attr(eval(gra), "gradient")
    Gras = as.matrix(apply(Gra, 2, sum))
    ss = matrix(rep(0, 25), nrow = 5, ncol = 5)
    hes = attr(eval(gra), "hessian")
    Hes = matrix(rep(0, 5 * 5), nrow = 5)
    for (i in 1:m) {
        Hes = Hes + hes[i, , ]
        ss = ss + Gra[i, ] \%*\% t(Gra[i, ])
    }
    vr = sqrt(diag(solve(ss)))
    im = -Hes
    IM = sqrt(1/m * diag(solve(im \%*\% t(im))))
    list(gd = Gras, hes = Hes, vr = vr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
