get_hessian <- function(y, theta, derivs)
{
    single_point_hessians <- lapply(y, function(yi){
        lambda <- theta$lambda
        mu1 <- theta$mu1
        mu2 <- theta$mu2
        var1 <- theta$var1
        var2 <- theta$var2
        drop(attr(eval(derivs),"hessian"))
    })
    hess <- Reduce("+", single_point_hessians)
    hess / length(y)
}

get_inverse_hessian <- function(y, theta, derivs)
{
    hess <- get_hessian(y, theta, derivs)
    solve(hess)
}

get_gradient <- function(y, theta, derivs)
{
    single_point_gradients <- lapply(y, function(yi){
        lambda <- theta$lambda
        mu1 <- theta$mu1
        mu2 <- theta$mu2
        var1 <- theta$var1
        var2 <- theta$var2
        drop(attr(eval(derivs),"grad"))
    })
    grad <- Reduce("+", single_point_gradients)
    grad / length(y)
}

mixture_newton <- function(y, maxit, tol, theta, verbose=TRUE)
{
    theta <- param(theta)
    
    log_lik_param <- expression(
        log(
            1/(1+exp(-lambda))*1/(2*pi*exp(var1))^(1/2)*
                exp(-(yi-mu1)*(yi-mu1)/(2*exp(var1))) +
                (1-1/(1+exp(-lambda)))*1/(2*pi*exp(var2))^(1/2)*
                exp(-(yi-mu2)*(yi-mu2)/(2*exp(var2)))    
        )
    )
    derivs <- deriv3(log_lik_param, c("lambda","mu1","mu2","var1","var2"))
    
    for(i in seq_len(maxit))
    {
        inv_hess <- tryCatch({
                get_inverse_hessian(y, theta, derivs)
            },
            error=function(cond) {
                message(paste("Ending Newton's search on iteration", i, 
                              "due to singular Hessian:"))
                message(cond)
                message("")
                return(NULL)
        })
        if (is.null(inv_hess))
        {
            break
        }
        grad <- get_gradient(y, theta, derivs)
        y_mat <- unlist(theta) - (inv_hess %*% grad)
        # Convert matrix with row names to named list
        new_theta <- as.list(data.frame(t(y_mat)))
        
        is_converged <- check_convergence(unparam(theta), unparam(new_theta), 
                                          tol=tol)
        
        theta <- new_theta
        
        if (is_converged)
        {
            break
        }
    }
    
    list(mle=get_return_theta(unparam(theta)), 
         stderr=get_stderr(y, unparam(theta)))
}