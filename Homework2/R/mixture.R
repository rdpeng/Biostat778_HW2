mixture <- function(y, method, maxit = NULL, tol = 1e-08, param0 = NULL) 
{
      y <- unlist(y)
      n <- length(y)
      method=match.arg(method,c("EM","newton"))
      
      ## Choose param0 if NULL is given
      ## parameter is of the form (lambda, mu1, mu2, sigma1, sigma2)
      ## note that parameters are estimated such that mu1 is always smaller than mu2
      
      if (is.null(param0))
      {
            ## find param0 in a naive way: by selecting peaks and
            ## calculate variance for sub-population around peak
            
            ## note that starting values picked do not guarantee
            ## convergence for newton method
            
            mu <- mean(y)
            mix <- density(y)
            if ((mix$x[mix$y == max(mix$y)])[1] <= mu)
            {
                  mixnorm1<-density(y[y <= mu])
                  mu10 <- mixnorm1$x[mixnorm1$y == max(mixnorm1$y)]
                  sigma10 <- mean((y[y <= mu10] - mu10) ^ 2)
                  
                  mixnorm2 <- density(y[y >= max(mu10 + sqrt(sigma10), mu)])
                  mu20 <- mixnorm2$x[mixnorm2$y == max(mixnorm2$y)]
                  sigma20 <- mean((y[y >= mu20] - mu20) ^ 2)
            }
            
            if ((mix$x[mix$y == max(mix$y)])[1] > mu)
            {
                  mixnorm2<-density(y[y >= mu])
                  mu20 <- mixnorm2$x[mixnorm2$y == max(mixnorm2$y)]
                  sigma20 <- mean((y[y >= mu20] - mu20) ^ 2)
                  
                  mixnorm1 <- density(y[y <= max(mu20 - sqrt(sigma20), mu)])
                  mu10 <- mixnorm1$x[mixnorm1$y == min(mixnorm1$y)]
                  sigma10 <- mean((y[y <= mu10] - mu10) ^ 2)
            }
            
            lambda0 <- (mu20 - mu) / (mu20 - mu10)
            param0 <- c(lambda0, mu10, mu20, sigma10, sigma20)
      }
      
      ## EM method for MLE
      EM <- function()
      {
            if (is.null(maxit))  maxit <- 500
            
            ## Calculate MLE
            for (i in 1:maxit)
            {
                  lambda0 <- param0[1]
                  mu10 <- param0[2]
                  mu20 <- param0[3]
                  sigma10 <- param0[4]
                  sigma20 <- param0[5]
                  
                  sd10 <- sqrt(sigma10)
                  sd20 <- sqrt(sigma20)
                  
                  theta <- lambda0 * dnorm(y, mu10, sd10) / 
                        (lambda0 * dnorm(y, mu10, sd10) + 
                               (1-lambda0) * dnorm(y, mu20, sd20))
                  
                  lambda <- sum(theta) / n
                  mu1 <- sum(theta * y) / sum(theta)
                  mu2 <- sum((1 - theta) * y) / sum(1 - theta)
                  sigma1 <- sum(theta * (y - mu1) ^ 2) / sum(theta)
                  sigma2 <- sum((1 - theta) * (y - mu2) ^ 2) / sum(1- theta)
                  
                  param <- c(lambda, mu1, mu2, sigma1, sigma2) 
                  
                  if (max(abs((param - param0) / param0)) > tol)
                        param0 <- param
                  else break
            }
            param
      }
      
      
      ## newton method for MLE
      newton <- function()
      {
            if (is.null(maxit))  maxit <- 100
            
            l <- expression(log(llambda * exp(-(y - lmu1) ^ 2 / (2 * lsigma1)) / sqrt(2 * pi * lsigma1) + 
                                      (1 - llambda) * exp(-(y - lmu2) ^ 2 / (2 * lsigma2)) / sqrt(2 * pi * lsigma2)))
            
            for (i in 1:maxit)
            {
                  lambda0 <- param0[1]
                  mu10 <- param0[2]
                  mu20 <- param0[3]
                  sigma10 <- param0[4]
                  sigma20 <- param0[5]
                  
                  der <- deriv3(l, c('llambda', 'lmu1', 'lmu2', 'lsigma1', 'lsigma2'))
                  llambda <- lambda0
                  lmu1 <- mu10
                  lmu2 <- mu20
                  lsigma1 <- sigma10
                  lsigma2 <- sigma20
                  
                  grad <- matrix(attr(eval(der), "gradient"), n, 5)
                  hess <- data.frame(attr(eval(der), "hessian"))
                  grad <- apply(grad, 2, sum)
                  hess <- matrix(apply(hess, 2, sum), 5, 5)
                  param <- as.vector(param0 - grad %*% solve(hess))
                  
                  if (max(abs((param - param0) / param0)) > tol)
                        param0 <- param
                  else break
            }
            param
      }
      
      ## Calculate mle for method
      if (method == "EM")     mle <- EM()
      if (method == "newton")       mle <- newton()
      
      ## Calculate fisher information
      l <- expression(log(llambda * exp(-(y - lmu1) ^ 2 / (2 * lsigma1)) / sqrt(2 * pi * lsigma1) + 
                                (1 - llambda) * exp(-(y - lmu2) ^ 2 / (2 * lsigma2)) / sqrt(2 * pi * lsigma2)))
      der <- deriv3(l, c('llambda', 'lmu1', 'lmu2', 'lsigma1', 'lsigma2'))
      llambda <- mle[1]
      lmu1 <- mle[2]
      lmu2 <- mle[3]
      lsigma1 <- mle[4]
      lsigma2 <- mle[5]
      
      grad <- attr(eval(der),"gradient")
      
      score <- matrix(0,5,5)
      for (i in 1:n)
      {
            score <- score + grad[i, ] %*% t(grad[i, ])
      }
      score <- score
      sd <- c(sqrt(diag(solve(score))))

      
      ## Return a list with elements `mle' for the maximum likelhood estimates and
      ## `stderr' for their standard errors.
      names(mle) <- c("lambda", "mu1", "mu2", "sigma1", "sigma2")
      names(sd) <- c("lambda", "mu1", "mu2", "sigma1", "sigma2")
      list(mle = mle, stderr = sd)
}