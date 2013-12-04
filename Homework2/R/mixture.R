mixture <- function(y, method=c("newton","EM"), maxit = NULL, tol = 1e-08, param0 = NULL) {
      method <- match.arg(method)
      if (is.null(param0))
            ###Newton's Method
            param0 <- c(0.5,mean(y),mean(y),sd(y),sd(y))
      if (method=="newton") {
            if (is.null(maxit))
                  maxit <- 100
            iteration <- 0
            tolerance <- tol + 1
            param <- param0
            while(iteration < maxit & tolerance >= tol) {
                  iteration <- iteration + 1
                  likederiv <- deriv3(~log(lambda*dnorm((y-mu1)/sigma1)/sigma1 + (1-lambda)*dnorm((y-mu2)/sigma2)/sigma2),c("lambda","mu1","mu2","sigma1","sigma2"),function(lambda,mu1,mu2,sigma1,sigma2,yi){})                                                    
                  derivres <- likederiv(param[1],param[2],param[3],param[4],param[5],y)
                  hessian <- colSums(attr(derivres,"hessian"))
                  oldparam <- param
                  param <- param - colSums(attr(derivres,"gradient")) %*% solve(hessian)
                  tolerance <- sum((oldparam-param)^2)
            }
            var <- diag(solve(-hessian)/length(y))
      } else if (method=="EM") {
            ###EM Algorithm
            if (is.null(maxit))
                  maxit <- 500
            iteration <- 0
            tolerance <- tol + 1
            param <- param0
            names(param) <- c("lambda", "mu1", "mu2", "sigma1", "sigma2")
            while(iteration < maxit & tolerance >= tol) {
                  iteration <- iteration + 1
                  oldparam <- param
                  bp <- oldparam[1]*dnorm(y,oldparam[2],oldparam[4])/(oldparam[1]*dnorm(y,oldparam[2],oldparam[4])+(1-oldparam[1])*dnorm(y,oldparam[3],oldparam[5]))
                  param[1] <- sum(bp)/length(bp)
                  param[2] <- sum(bp*y)/sum(bp)
                  param[3] <- sum((1-bp)*y)/sum((1-bp))
                  param[4] <- sqrt(sum(bp*(y-param[2])^2)/sum(bp))
                  param[5] <- sqrt(sum((1-bp)*(y-param[3])^2)/sum(1-bp))
                  tolerance <- sum((oldparam-param)^2)
            }        
            fisherinf <- matrix(0,nrow=5,ncol=5)
            for (i in 1:length(y)) {
                  score <- c( (bp[i]-param[1])/(param[1]*(1-param[1])),                ##lambda
                              bp[i]*(y[i]-param[2])/(param[4]^2),                      ##mu1       
                              (1-bp[i])*(y[i]-param[3])/(param[5]^2),                  ##mu2
                              bp[i]*(-1/param[4]+(y[i]-param[2])^2/param[4]^3),        ##sigma1
                              (1-bp[i])*(-1/param[5]+(y[i]-param[3])^2/param[5]^3)     ##sigma1
                  )
                  fisherinf <- fisherinf + score %*% t(score)
            }
            var <- diag(solve(fisherinf/length(y)))/length(y)
      } else {
            stop("Please give a valid method.")
      }
      names(var) <- c("lambda", "mu1", "mu2", "sigma1", "sigma2")
      return(list(mle=param,stderr=sqrt(var)))    
}

