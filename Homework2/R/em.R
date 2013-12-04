get_expected_z <- function(y, theta)
{
    lambda <- theta$lambda
    mu1 <- theta$mu1
    mu2 <- theta$mu2
    sigma1 <- sqrt(theta$var1)
    sigma2 <- sqrt(theta$var2)
    # BAYES' rule
    # P(N1 | Y=y) = P(Y=y | N1)*P(N1) / (P(Y=y | N1)*P(N1) + P(Y=y | N2)*P(N2))
    dnorm(y, mean=mu1, sd=sigma1)*lambda/ 
        (lambda*dnorm(y, mean=mu1, sd=sigma1) + 
            (1-lambda)*dnorm(y, mean=mu2, sd=sigma2))
}

get_theta <- function(y, z)
{
    theta <- list()
    theta$lambda <- mean(z)
    theta$mu1 <- sum(z*y)/sum(z)
    theta$mu2 <- sum((1-z)*y)/sum(1-z)
    theta$var1 <- sum(z*(y - theta$mu1)*(y - theta$mu1))/sum(z)
    theta$var2 <- sum((1-z)*(y - theta$mu2)*(y - theta$mu2))/sum(1-z)
    theta
}

mixture_em <- function(y, maxit, tol, theta, verbose=TRUE)
{
    z <- get_expected_z(y, theta)
    
    for(i in seq_len(maxit))
    {
        new_theta <- get_theta(y, z)
        z <- get_expected_z(y, new_theta)
        
        is_converged <- check_convergence(theta, new_theta, tol=tol)
        
        theta <- new_theta
        
        if (is_converged)
        {
            break
        }
    }
    
    list(mle=get_return_theta(theta), stderr=get_stderr(y, theta))
}
