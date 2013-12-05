pkgname <- "Homework2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Homework2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Homework2-package")
### * Homework2-package

flush(stderr()); flush(stdout())

### Name: Homework2-package
### Title: Homework 2 for Advanced Computing
### Aliases: Homework2-package Homework2

### ** Examples

w = 0.6 
m1 = 5; m2 = -5; v1 = 2; v2 = 3

W = rbinom(1000,1,prob=w)
N1 = rnorm(1000,m1,sqrt(v1)) 
N2 = rnorm(1000,m2,sqrt(v2))
Y = W*N1+(1-W)*N2
rm(w,m1,m2,v1,v2,W,N1,N2)

mixture(y=Y,method="Newton",maxit = 100)
mixture(y=Y,method="EM",maxit = 300)



cleanEx()
nameEx("mixture")
### * mixture

flush(stderr()); flush(stdout())

### Name: mixture
### Title: Estimating Gaussian mixture model parameters
### Aliases: mixture

### ** Examples

w = 0.6 
m1 = 5; m2 = -5; v1 = 2; v2 = 3

W = rbinom(1000,1,prob=w)
N1 = rnorm(1000,m1,sqrt(v1)) 
N2 = rnorm(1000,m2,sqrt(v2))
Y = W*N1+(1-W)*N2
rm(w,m1,m2,v1,v2,W,N1,N2)

mixture(y=Y,method="Newton",maxit = 100)
mixture(y=Y,method="EM",maxit = 300)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
