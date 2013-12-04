get_hessian <- function(y, theta, derivs)
{
    single_point_hessians <- lapply(y, function(yi){
        lambda <- theta$lambda
        mu1 <- theta$mu1
        mu2 <- theta$mu2
        sigma1 <- theta$sigma1
        sigma2 <- theta$sigma2
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
        sigma1 <- theta$sigma1
        sigma2 <- theta$sigma2
        drop(attr(eval(derivs),"grad"))
    })
    grad <- Reduce("+", single_point_gradients)
    grad / length(y)
}

get_hess_stderr <- function(y, theta)
{
    derivs <- get_derivs()
    inv_hess <- tryCatch({
        inv_hess <- get_inverse_hessian(y, theta, derivs)
        inv_hess_det <- det(inv_hess)
        if (abs(inv_hess_det) > 1e12)
        {
            stop(paste("Hessian is computationally singular: det =",
                       1/inv_hess_det))
        }
        inv_hess
    },
    error=function(cond) {
        message(cond)
        message()
        return(NULL)
    })
    if (is.null(inv_hess))
    {
        return(get_error_vec())
    } else 
    {
        stderr <- sqrt(-diag(inv_hess))
        return(stderr)
    }
}

mixture_newton <- function(y, maxit, tol, theta, verbose=TRUE)
{
    error <- FALSE
    
    derivs <- get_derivs_param()
    theta <- param(theta)
    
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
            error <- TRUE
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
    
    mle <- get_return_theta(unparam(theta))
    stderr <- get_error_vec()
    if (!error)
    {
        stderr <- get_hess_stderr(y, unparam(theta))
    }
    
    list(mle=mle, 
         stderr=stderr)
}