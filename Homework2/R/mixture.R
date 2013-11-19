mixture <- function(y, method = c("newton", "EM"), maxit = NULL, tol = 1e-8,
                    param0 = NULL) {
        method <- match.arg(method)
        switch(method,
               newton = {
                       maxit <- 100
                       newton.mix(y, param0, maxit, tol)
               },
               EM = {
                       maxit <- 500
                       em.mix(y, param0, maxit)
               })
}

newton.mix <- function(y, param0, maxit, tol) {
        grad <- deriv3(~ log(lambda * dnorm((x - m1)/s1)/s1 + (1-lambda) * dnorm((x-m2)/s2)/s2), c("lambda", "m1", "m2", "s1", "s2"), c("y", "lambda", "m1", "m2", "s1", "s2"))
        drv <- function(y, param) {
                g <- grad(y, param[1], param[2], param[3], param[4], param[5])
                list(likelihood = sum(g),
                     gradient = colSums(attr(g, "gradient")),
                     hessian = colSums(attr(g, "hessian"), dims = 1))
        }
        d <- drv(y, param0)
        ll0 <- d$likelihood
        grad <- d$gradient
        hess <- d$hessian
        for(i in seq_len(maxit)) {
                param <- solve(hess, hess %*% param0 - grad)
                d <- drv(y, param)
                ll <- d$likelihood
                delta <- abs(ll - ll0)
                if(delta < tol) {
                        convergence <- 0
                        break
                }
                else {
                        ll0 <- ll
                        grad <- d$gradient
                        hess <- d$hessian
                }
        }
        if(i == maxit && delta >= tol)
                convergence <- 1
        hess <- drv(y, param)$hessian
        list(mle = param, stderr = sqrt(diag(hess)), convergence = convergence)
}

em.mix <- function(y, param0, maxit, tol) {

}
