# Useful functions

get_lik_one <- function(theta, yi)
{
    lambda <- theta[[1]]
    mu1 <- theta[[2]]
    mu2 <- theta[[3]]
    sigma1 <- theta[[4]]
    sigma2 <- theta[[5]]
    log(
        lambda/(2*pi*sigma1^2)^(1/2)*exp(-(yi-mu1)*(yi-mu1)/(2*sigma1^2)) +
        (1-lambda)/(2*pi*sigma2^2)^(1/2)*exp(-(yi-mu2)*(yi-mu2)/(2*sigma2^2))    
    )
}

get_lik <- function(theta, y)
{
    lik_ones <- sapply(y, get_lik_one, theta=theta)
    sum(lik_ones)
}

get_initial_theta <- function(y)
{
    d <- density(y)
    maxima <- which(diff(sign(diff(d$y))) == -2)
    arg_max_x <- d$x[maxima+1]
    
    if (length(arg_max_x) > 2)
    {
        largest <- order(-d$y[maxima+1])[1:2]
        arg_max_x <- arg_max_x[largest]
    }
    
    theta <- list()
    
    if (length(arg_max_x) == 2)
    {
        # Get the minimum between two maxima
        minima <- which(diff(sign(diff(d$y))) == 2)
        arg_min_x <- d$x[minima+1]
        middle_mins <- arg_min_x[arg_min_x > arg_max_x[1] & 
                                 arg_min_x < arg_max_x[2]]
        # Dividing line in between maxima
        div_min <- middle_mins[1]
        
        # Density 1
        f1 <- y[y <= div_min]
        # Density 2
        f2 <- y[y > div_min]
        theta$lambda <- length(f1)/length(y)
        theta$mu1 <- arg_max_x[1]
        theta$mu2 <- arg_max_x[2]
        theta$sigma1 <- sd(f1)
        theta$sigma2 <- sd(f2)
    } else if (length(arg_max_x) == 1)
    {
        theta$lambda <- 0.5
        theta$mu1 <- arg_max_x - sqrt(var(y))*0.05
        theta$mu2 <- arg_max_x + sqrt(var(y))*0.05
        theta$sigma1 <- sd(y)
        theta$sigma2 <- sd(y)
    } else 
    {
        stop("Density of y has no maximum. Strange things are happening!")
    }
    theta
}

check_convergence <- function(old_theta, new_theta, tol=1e-8, 
                              method=c("l_2","l_inf"))
{
    method <- match.arg(method)
    delta <- sapply(1:length(old_theta), function(x) {
        abs(old_theta[[x]] - new_theta[[x]])
    })
    if (method == "l_inf")
    {
        return(all(delta < tol))
    } else if (method == "l_2")
    {
        return(sum(delta*delta) < tol*tol)
    }
    else {
        stop(paste("unknown distance metric:",method))
    }
}

get_derivs <- function()
{
    # Observed data likelihood
    log_lik <- expression(
        log(
            lambda/(2*pi*sigma1^2)^(1/2)*
            exp(-(yi-mu1)*(yi-mu1)/(2*sigma1^2)) +
            (1-lambda)/(2*pi*sigma2^2)^(1/2)*
            exp(-(yi-mu2)*(yi-mu2)/(2*sigma2^2))    
        )
    )
    derivs <- deriv3(log_lik, c("lambda","mu1","mu2","sigma1","sigma2"))
    return(derivs)
}

get_derivs_param <- function()
{
    log_lik_param <- expression(
        log(
            1/(1+exp(-lambda))*1/(2*pi*exp(2*sigma1))^(1/2)*
                exp(-(yi-mu1)*(yi-mu1)/(2*exp(2*sigma1))) +
                (1-1/(1+exp(-lambda)))*1/(2*pi*exp(2*sigma2))^(1/2)*
                exp(-(yi-mu2)*(yi-mu2)/(2*exp(2*sigma2)))    
        )
    )
    derivs <- deriv3(log_lik_param, c("lambda","mu1","mu2","sigma1","sigma2"))
    return(derivs)
}

get_information <- function(y, theta)
{
    # Observed data likelihood
    derivs <- get_derivs()
    
    single_point_scores <- lapply(y, function(yi){
        lambda <- theta$lambda
        mu1 <- theta$mu1
        mu2 <- theta$mu2
        sigma1 <- theta$sigma1
        sigma2 <- theta$sigma2
        score <- drop(attr(eval(derivs),"grad"))
        tcrossprod(score)
    })
    score_sum <- Reduce("+", single_point_scores)
    score_sum / length(y)
}

get_error_vec <- function()
{
    vec <- rep(NaN, 5)
    names(vec) <- c("lambda", "mu1", "mu2", "sigma1", "sigma2")
    return(vec)
}

get_stderr <- function(y, theta)
{
    I <- get_information(y, theta)
    I_inv <- tryCatch(
        solve(I), 
        error=function(cond){
            message("Information matrix is singular")
            return(NULL)
    })
    if (is.null(I_inv))
    {
        return(get_error_vec())
    } else 
    {
        stderr <- diag(I_inv)
        names(stderr) <- c("lambda", "mu1", "mu2", "sigma1", "sigma2")
        return(sqrt(stderr))
    }
}

get_theta_from_vec <- function(theta)
{
    names(theta) <- NULL
    new_theta <- list(lambda=theta[1], mu1=theta[2], mu2=theta[3], 
                      sigma1=theta[4], sigma2=theta[5])
    new_theta
}

get_return_theta <- function(theta)
{
    new_theta <- unlist(theta)
    names(new_theta) <- c("lambda","mu1","mu2","sigma1","sigma2")
    new_theta
}

param <- function(theta)
{
    theta$lambda <- log(theta$lambda/(1-theta$lambda))
    theta$sigma1 <- log(theta$sigma1)
    theta$sigma2 <- log(theta$sigma2)
    theta
}

unparam <- function(theta)
{
    theta$lambda <- 1/(1+exp(-theta$lambda))
    theta$sigma1 <- exp(theta$sigma1)
    theta$sigma2 <- exp(theta$sigma2)
    theta
}