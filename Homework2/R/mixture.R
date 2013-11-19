mixture <- function(y, method = c("newton", "EM"), param0 = NULL) {
        method <- match.arg(method)
        switch(method,
               newton = newton.mix(y, param0),
               EM = em.mix(y, param0))
}

newton.mix <- function(y, param0) {
        
}

em.mix <- function(y, param0) {

}
