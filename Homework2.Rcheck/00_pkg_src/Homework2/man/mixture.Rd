\name{mixture}
\alias{mixture}
\alias{EM.Newton.Mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
        EM and Newton for Mixture Normal Parameters Estimation
%%  ~~function to do ... ~~
}
\description{
        Use EM or Newton algorithm to estimate parameters in the mixture normal model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mixture(y, method, maxit = NULL, tol = 1e-08, param0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
        y is a vector of the observed data from a mixture normal model.
%%     ~~Describe \code{y} here~~
}
  \item{method}{
        method can take EM or newton depending on which algorithm people want to use.
%%     ~~Describe \code{method} here~~
}
  \item{maxit}{
        maxit is the max number of iteration and default value is 100 for Newton and 500 for EM.
%%     ~~Describe \code{maxit} here~~
}
  \item{tol}{
        tol is the tolerance to control the convergence.
%%     ~~Describe \code{tol} here~~
}
  \item{param0}{
        param0 is the initial parameters input for the model and in case people didn't enter any initial value, there is default specified in the program.
%%     ~~Describe \code{param0} here~~
}
}
\details{
        Use EM or Newton algorithm to estimate parameters in the mixture normal model.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
        \item{mle}{
                return the mle value for the five parameters in the mixed two normal model.
        }
        \item{stderr}{
                return the asymptotic standard error for the parameters in the mixed two normal model.
        }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
        Advanced Statistics Computing Class, Dr. Peng.
%% ~put references to the literature/web site here ~
}
\author{
        Yu Du
%%  ~~who you are~~
}
\examples{
z=numeric(0)
y=numeric(0)
for(i in 1:1000){
        z[i]=rbinom(1,1,0.2)
        if(z[i]==1){
                y[i]=rnorm(1,10,1)
        }else{
                y[i]=rnorm(1,3,1)
        }
}
mixture(y,"EM")


function (y, method, maxit = NULL, tol = 1e-08, param0 = NULL) 
{
    y = unlist(y)
    n = length(y)
    method = match.arg(method, c("EM", "newton"))
    if (method == "EM") {
        e = numeric(0)
        tlam = numeric(0)
        tmu1 = numeric(0)
        tmu2 = numeric(0)
        tsigma1 = numeric(0)
        tsigma2 = numeric(0)
        if (is.null(param0)) {
            tlam[1] = 0.1
            tmu1[1] = 3
            tmu2[1] = 2
            tsigma1[1] = 4
            tsigma2[1] = 5
        }
        else {
            tlam[1] = param0[1]
            tmu1[1] = param0[2]
            tmu2[1] = param0[3]
            tsigma1[1] = param0[4]
            tsigma2[1] = param0[5]
        }
        if (!is.numeric(maxit)) {
            maxit = 500
        }
        for (j in 1:maxit) {
            e = tlam[j] * dnorm(y, tmu1[j], tsigma1[j]^(1/2))/(tlam[j] * 
                dnorm(y, tmu1[j], tsigma1[j]^(1/2)) + (1 - tlam[j]) * 
                dnorm(y, tmu2[j], tsigma2[j]^(1/2)))
            tlam[j + 1] = mean(e)
            tmu1[j + 1] = sum(y * e)/sum(e)
            tmu2[j + 1] = sum(y * (1 - e))/sum(1 - e)
            tsigma1[j + 1] = sum(e * ((y - tmu1[j + 1])^2))/sum(e)
            tsigma2[j + 1] = sum((1 - e) * ((y - tmu2[j + 1])^2))/sum(1 - 
                e)
            if (sqrt((tlam[j + 1] - tlam[j])^2 + (tmu1[j + 1] - 
                tmu1[j])^2 + (tmu2[j + 1] - tmu2[j])^2 + (tsigma1[j + 
                1] - tsigma1[j])^2 + (tsigma2[j + 1] - tsigma2[j])^2) <= 
                tol) {
                break
            }
        }
        lambda = tlam[j + 1]
        mu1 = tmu1[j + 1]
        mu2 = tmu2[j + 1]
        sigma1 = tsigma1[j + 1]
        sigma2 = tsigma2[j + 1]
        l = expression(log(lambda * exp(-(y - mu1)^2/(2 * sigma1))/sqrt(2 * 
            pi * sigma1) + (1 - lambda) * exp(-(y - mu2)^2/(2 * 
            sigma2))/sqrt(2 * pi * sigma2)))
        de = deriv3(l, c("lambda", "mu1", "mu2", "sigma1", "sigma2"))
        lambda = lambda
        mu1 = mu1
        mu2 = mu2
        sigma1 = sigma1
        sigma2 = sigma2
        grad = attr(eval(de), "gradient")
        score = matrix(0, 5, 5)
        for (i in 1:n) {
            score = score + grad[i, ] \%*\% t(grad[i, ])
        }
        cov = solve(score)
        se = c(sqrt(diag(cov)))
        return(list(mle = c(lambda = lambda, mu1 = mu1, mu2 = mu2, 
            sigma1 = sigma1, sigma2 = sigma2), stderr = c(lamda = se[1], 
            mu1 = se[2], mu2 = se[3], sigma1 = se[4], sigma2 = se[5])))
    }
    if (method == "newton") {
        if (is.null(param0)) {
            temp = mixture(y, "EM")[[1]]
            add = runif(4, -0.5, 0.5)
            add = as.vector(cbind(runif(1, -0.05, 0.05), t(add)))
            theta = temp + add
        }
        else {
            theta = param0
        }
        if (!is.numeric(maxit)) {
            maxit = 100
        }
        l = expression(log(lambda * exp(-(y - mu1)^2/(2 * sigma1))/sqrt(2 * 
            pi * sigma1) + (1 - lambda) * exp(-(y - mu2)^2/(2 * 
            sigma2))/sqrt(2 * pi * sigma2)))
        de = deriv3(l, c("lambda", "mu1", "mu2", "sigma1", "sigma2"))
        for (i in 1:maxit) {
            lambda = theta[1]
            mu1 = theta[2]
            mu2 = theta[3]
            sigma1 = theta[4]
            sigma2 = theta[5]
            grad = attr(eval(de), "gradient")
            grad = matrix(apply(grad, 2, sum), 5, byrow = T)
            hes = attr(eval(de), "hessian")
            hes = data.frame(hes)
            hes = matrix(apply(hes, 2, sum), 5, byrow = T)
            theta = theta - solve(hes) \%*\% grad
            if (sqrt(sum((solve(hes) \%*\% grad)^2)) <= tol) {
                break
            }
        }
        lambda = theta[1]
        mu1 = theta[2]
        mu2 = theta[3]
        sigma1 = theta[4]
        sigma2 = theta[5]
        grad = attr(eval(de), "gradient")
        score = matrix(0, 5, 5)
        for (i in 1:n) {
            score = score + grad[i, ] \%*\% t(grad[i, ])
        }
        cov = solve(score)
        se = c(sqrt(diag(cov)))
        return(list(mle = c(lamda = theta[1], mu1 = theta[2], 
            mu2 = theta[3], sigma1 = theta[4], sigma2 = theta[5]), 
            stderr = c(lamda = se[1], mu1 = se[2], mu2 = se[3], 
                sigma1 = se[4], sigma2 = se[5])))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
