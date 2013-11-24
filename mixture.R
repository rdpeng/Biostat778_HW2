### Estimate the parameters from a 2 component guassian mixture model.
### Options for using either EM algorith or Newton's method.

mixture <- function(y, method, maxit=NULL, tol = 1e-08, param0 = NULL) {

    ## If initial values not provided, estimate using k-means clustering.
    if(is.null(param0)) {
        v <- kmeans(y, 2, nstart=10)
        mu0 <- sort(v$centers)
        s20 <- v$withinss/(v$size - 1)[order(v$centers)]
        nn <- pars$size[order(pars$centers)]
        lambda0 <- nn/sum(nn)
    }

    ## Newton's method ##

    ## EM algorithm ##
    mu <- mu0
    s2 <- s20
    lambda <- lambda0
    ## use repeat{...} instead of while?
    while(1) {
        ## E STEP ##
        z.tilde <- lambda[1]*dnorm(y, mu[1], s2[1]) /
        (lambda[1]*dnorm(y, mu[1], s2[1]) + lambda[2]*dnorm(y, mu[2], s2[2]))
            ## return a list with elements `mle' and `stderr'.

        ## M STEP ##
        mu[1] <- sum(z.tilde * y) / sum(z.tilde)
        mu[2] <- sum((1-z.tilde) * y) / sum(1-z.tilde)
        s2[1] <- sum(z.tilde * (y-mu[1])^2 ) / sum(z.tilde)
        s2[2] <- sum((1-z.tilde) * (y-mu[2])^2 ) / sum(1-z.tilde)
        lambda <- sum(z.tilde)/length(z.tilde)

        if(something < tol) return(0) # terminate loop (crossprod of estimates?)
    }
    ## find standard errors

}


## generate test data
library(mixtools)
set.seed(1234)
n <- 1000
pi <- c(.3, .7)
mu <- c(0, 2.5)
s2 <- c(.3, 1)
y <- rnormmix(n, pi, mu, sqrt(s2))
hist(y, breaks=50, col='gray', border='gray')
