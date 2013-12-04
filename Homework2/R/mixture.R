#' Mixture density determination
#'
#' @param y a numeric vector containing the observations from the mixture 
#' density
#' @param method either \code{"EM"} for EM algorithm or \code{"newton"} for 
#' Newton's method
#' @param maxit the maximum number of iterations for the method to run. If 
#' \code{NULL} then the method defaults to 500 for EM and 100 for Newton's 
#' method.
#' @param tol the tolerance to determine convergence: if the method updates the 
#' parameter and the new value is less than \code{tol} away in Euclidean 
#' distance, then convergence is achieved.
#' @param param0 a numeric vector, the starting value for the parameter. Should 
#' be in the following order: \code{lambda}, \code{mu1}, \code{mu2}, 
#' \code{sigma1}, \code{sigma2}. If \code{NULL}, then the method attempts to 
#' determine initial values.
#' @return \code{mixture()} returns a list with the following components: 
#' \item{mle}{a numeric vector containing the last update of the parameter.} 
#' \item{stderr}{a numeric vector containing the standard errors of the 
#' parameters from the observed information.}
#' @export
#' @examples
#' set.seed(2013-12-2)
#' y <- c(rnorm(100,mean=-2), rnorm(100, mean=2))
#' theta <- list(lambda=0.5, mu1=-2, mu2=2, simga1=1, sigma2=1)
#' mixture(y, "newton", param0=unlist(theta))
#' mixture(y, "EM", param0=unlist(theta))
#' \dontrun{
#' head(hw2_data, 20)
#' mixture(hw2_data, "EM")
#' }
mixture <- function(y, method, maxit = NULL, tol = 1e-08, param0 = NULL) {
    if (is.null(param0))
    {
        param0 <- get_initial_theta(y)
    } else {
        param0 <- get_theta_from_vec(param0)
    }
    
    if (method == "newton")
    {
        if (is.null(maxit))
        {
            maxit <- 100
        }
        return(mixture_newton(y, maxit, tol, param0))
    } else if (method == "EM")
    {
        if (is.null(maxit))
        {
            maxit <- 500
        }
        return(mixture_em(y, maxit, tol, param0))
    } else {
        stop(paste("Unrecognized method:",method))
    }
}