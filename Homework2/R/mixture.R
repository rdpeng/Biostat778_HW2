mixture = function(y, method, maxit = NULL, tol = 1e-08, param0 = NULL) {
     y = unlist(y)
     n = length(y)
     iters = 0
     diff = tol+5
     
     tryCatch( {
          method = match.arg(method,choices=c("newton","EM"))
     }, error = function(err) {
          stop("Invalid method choice: choose either newton or EM")
     })
     
     if (is.null(maxit)) {
          if (method == "newton")
               maxit = 100
          if (method == "EM")
               maxit = 500
     }
     if (is.null(param0)) {
          # Choose initial parameters based on k-means clustering
          # NOTE: Sometimes clusters 1 and 2 are switched by kmeans
          # So mu1, mu2 and sigma1,sigma2 (and their std errors)
          # switch values when running this function at different times
          # This initialization may not be numerically stable
          clusters = kmeans(y,centers=2)
          lambda = clusters$size[1]/n
          mu1 = clusters$center[1]
          mu2 = clusters$center[2]
          clust_sd = sqrt(clusters$withinss/clusters$size)
          sigma1 = clust_sd[1]
          sigma2 = clust_sd[2]
          theta = matrix(c(lambda,mu1,mu2,sigma1,sigma2),nrow=5,ncol=1)
     }
     
     # Minimize negative log-likelihood 
     if (method == "newton") {
          loglikterm = expression(-log((lambda/sigma1)*exp(-(y-mu1)^2/(sigma1^2)) 
                                       + ((1-lambda)/sigma2)*exp(-(y-mu2)^2/(sigma2^2))))
          dl = deriv3(loglikterm,namevec=c("lambda","mu1","mu2","sigma1","sigma2"))
          while (diff > tol && iters < maxit) {
               iters = iters + 1
               vals_dl = eval(dl)
               tryCatch({
                    grad = attr(vals_dl,"gradient")
                    grad = colSums(grad)
               }, error = function(err) {
                    stop("Unable to compute gradient")
               })
#                gradlam = sum(grad[,1])
#                gradmu1 = sum(grad[,2])
#                gradmu2 = sum(grad[,3])
#                gradsig1 = sum(grad[,4])
#                gradsig2 = sum(grad[,5])
               tryCatch({
               hess = attr(vals_dl,"hessian")
               hess = colSums(hess)
               }, error = function(err) {
                    stop("Unable to compute hessian")
               })
#                hesslam = sum(hess[1:n])
#                hessmu1 = sum(hess[(6*n+1):(7*n)])
#                hessmu2 = sum(hess[(12*n+1):(13*n)])
#                hesssig1 = sum(hess[(18*n+1):(19*n)])
#                hesssig2 = sum(hess[(24*n+1):(25*n)])
               
               # Newton directions
               tryCatch({
                    p = solve(hess)%*%grad
               }, error = function(err) {
                    stop("Could not compute inverse of Hessian")
               })
#                plam = -(1/hesslam)*gradlam
#                pmu1 = -(1/hessmu1)*gradmu1
#                pmu2 = -(1/hessmu2)*gradmu2
#                psigma1 = -(1/hesssig1)*gradsig1
#                psigma2 = -(1/hesssig2)*gradsig2
               
               # Find new theta
               theta_new = theta + p
#                lambda = lambda + plam
#                mu1 = mu1 + pmu1
#                mu2 = mu2 + pmu2
#                sigma1 = sigma1 + psigma1
#                sigma2 = sigma2 + psigma2
               
               # Calculate diff and update
#                diff = plam + pmu1 + pmu2 + psigma1 + psigma2
               diff = sum(abs(theta_new-theta))
               theta = theta_new
               lambda = theta[1]
               mu1 = theta[2]
               mu2 = theta[3]
               sigma1 = theta[4]
               sigma2 = theta[5]
               
               # Stop if NaNs are produced
               if (is.nan(diff)) {
                    mle = c(lambda,mu1,mu2,sigma1,sigma2)
                    names(mle) = c("lambda","mu1","mu2","sigma1","sigma2")
                    stderr = c(NA,NA,NA,NA,NA)
                    names(stderr) = c("lambda","mu1","mu2","sigma1","sigma2")
                    return(list(mle=mle,stderr=stderr))
                    stop("Error: NaN's were produced")
               }
          } # end while
          # Calculate standard errors
          # Asymptotic variance of theta = inverse of Fisher information matrix
          vals_dl = eval(dl)
          hess = attr(vals_dl,"hessian")
          hess = colSums(hess) # No negative here because negative hessian for negative log-lik-->positive hessian
          theta_var = solve(hess)
     }
     
     if (method == "EM") {
          # EM algorithm
          while (diff > tol && iters < maxit) {
               iters = iters + 1
               # E-step: Calculate Q function
               norm1 = dnorm(y,mu1,sigma1)
               norm1log = dnorm(y,mu1,sigma1,log=TRUE)
               norm2 = dnorm(y,mu2,sigma2)
               norm2log = dnorm(y,mu2,sigma2,log=TRUE)
               condexp_zi = lambda*norm1/(lambda*norm1+(1-lambda)*norm2)
               sumcondexp_zi = sum(condexp_zi)
               Q = sum(condexp_zi*norm1log + (1-condexp_zi)*norm2log)
               # M-step: Maximize Q function
               lambda_new = mean(condexp_zi)
               mu1_new = sum(condexp_zi*y)/sumcondexp_zi
               mu2_new = sum((1-condexp_zi)*y)/(n-sumcondexp_zi)
               sigma1_new = sqrt(sum(condexp_zi*((y-mu1_new)^2))/sumcondexp_zi)
               sigma2_new = sqrt(sum((1-condexp_zi)*((y-mu2_new)^2))/(n-sumcondexp_zi))
               # Calculate difference and update parameters
               diff = lambda_new-lambda + mu1_new-mu1 + mu2_new-mu2 + sigma1_new-sigma1 + sigma2_new-sigma2
               lambda = lambda_new
               mu1 = mu1_new
               mu2 = mu2_new
               sigma1 = sigma1_new
               sigma2 = sigma2_new
          } # end while
          # Calculate standard errors
          loglikterm = expression(-log((lambda/sigma1)*exp(-(y-mu1)^2/(sigma1^2)) 
                                       + ((1-lambda)/sigma2)*exp(-(y-mu2)^2/(sigma2^2))))
          dl = deriv3(loglikterm,namevec=c("lambda","mu1","mu2","sigma1","sigma2"))
          vals_dl = eval(dl)
          grad = attr(vals_dl,"gradient")
          fisherinfo = t(grad)%*%grad/n
          theta_var = solve(fisherinfo)
     }
     mle = c(lambda,mu1,mu2,sigma2,sigma2)
     names(mle) = c("lambda","mu1","mu2","sigma1","sigma2")
     stderr = sqrt(diag(theta_var))
     names(stderr) = c("lambda","mu1","mu2","sigma1","sigma2")
     
     return(list(mle=mle,stderr=stderr))
}