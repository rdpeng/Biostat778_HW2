# Useful functions

get_lik_one <- function(theta, yi)
{
    lambda <- theta[[1]]
    mu1 <- theta[[2]]
    mu2 <- theta[[3]]
    var1 <- theta[[4]]
    var2 <- theta[[5]]
    log(
        lambda/(2*pi*var1)^(1/2)*exp(-(yi-mu1)*(yi-mu1)/(2*var1)) +
        (1-lambda)/(2*pi*var2)^(1/2)*exp(-(yi-mu2)*(yi-mu2)/(2*var2))    
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
        theta$var1 <- var(f1)
        theta$var2 <- var(f2)
    } else if (length(arg_max_x) == 1)
    {
        theta$lambda <- 0.5
        theta$mu1 <- arg_max_x - sqrt(var(y))*0.05
        theta$mu2 <- arg_max_x + sqrt(var(y))*0.05
        theta$var1 <- var(y)
        theta$var2 <- var(y)
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

get_information <- function(y, theta)
{
    # Observed data likelihood
    log_lik <- expression(
        log(
            lambda/(2*pi*s1^2)^(1/2)*exp(-(yi-mu1)*(yi-mu1)/(2*s1^2)) +
                (1-lambda)/(2*pi*s2^2)^(1/2)*exp(-(yi-mu2)*(yi-mu2)/(2*s2^2))    
        )
    )
    derivs <- deriv3(log_lik, c("lambda","mu1","mu2","s1","s2"))
    
    single_point_scores <- lapply(y, function(yi){
        lambda <- theta$lambda
        mu1 <- theta$mu1
        mu2 <- theta$mu2
        s1 <- sqrt(theta$var1)
        s2 <- sqrt(theta$var2)
        score <- drop(attr(eval(derivs),"grad"))
        tcrossprod(score)
    })
    score_sum <- Reduce("+", single_point_scores)
    score_sum / length(y)
}

get_stderr <- function(y, theta)
{
    I <- get_information(y, theta)
    I_inv <- tryCatch(
        solve(I), 
        error=function(cond){
            stop("Information matrix is singular")
    })
    stderr <- diag(I_inv)
    names(stderr) <- c("lambda", "mu1", "mu2", "sigma1", "sigma2")
    sqrt(stderr)
}

get_theta_from_vec <- function(theta)
{
    names(theta) <- NULL
    new_theta <- list(lambda=theta[1], mu1=theta[2], mu2=theta[3], 
                      var1=theta[4]*theta[4], var2=theta[5]*theta[5])
    new_theta
}

get_return_theta <- function(theta)
{
    new_theta <- unlist(theta)
    new_theta[4] <- sqrt(new_theta[4])
    new_theta[5] <- sqrt(new_theta[5])
    names(new_theta) <- c("lambda","mu1","mu2","sigma1","sigma2")
    new_theta
}

param <- function(theta)
{
    theta$lambda <- log(theta$lambda/(1-theta$lambda))
    theta$var1 <- log(theta$var1)
    theta$var2 <- log(theta$var2)
    theta
}

unparam <- function(theta)
{
    theta$lambda <- 1/(1+exp(-theta$lambda))
    theta$var1 <- exp(theta$var1)
    theta$var2 <- exp(theta$var2)
    theta
}