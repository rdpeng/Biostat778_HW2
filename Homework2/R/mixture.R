mixture = function(y, method, maxit = NULL, tol = 1e-08, param0 = NULL) {
     y = unlist(y)
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
          # Choose initial parameters
          clusters = kmeans(y,centers=2)
          mu1 = clusters$center[1]
          mu2 = clusters$center[2]
          clust_sd = sqrt(clusters$withinss/clusters$size)
          sigma1 = clust_sd[1]
          sigma2 = clust_sd[2]
          theta = matrix(c(lambda,mu1,mu2,sigma1,sigma2),nrow=5,ncol=1)
     }
     n = length(y)
     iters = 0
     diff = tol+5
     
     # Minimize negative log-likelihood 
     if (method == "newton") {
          mu1 = 5
          mu2 = 10
          sigma1 = 1
          sigma2 = 2
          # Newton's method
          loglikterm = expression(-log((lambda/sigma1)*exp(-(y-mu1)^2/(sigma1^2)) 
                                       + ((1-lambda)/sigma2)*exp(-(y-mu2)^2/(sigma2^2))))
          dl = deriv3(loglikterm,namevec=c("lambda","mu1","mu2","sigma1","sigma2"))
          while (diff > tol && iters < maxit) {
               iters = iters + 1
               vals_dl = eval(dl)
               grad = attr(vals_dl,"gradient")
               grad = colSums(grad)
#                gradlam = sum(grad[,1])
#                gradmu1 = sum(grad[,2])
#                gradmu2 = sum(grad[,3])
#                gradsig1 = sum(grad[,4])
#                gradsig2 = sum(grad[,5])
               hess = attr(vals_dl,"hessian")
               hess = colSums(hess)
#                hesslam = sum(hess[1:n])
#                hessmu1 = sum(hess[(6*n+1):(7*n)])
#                hessmu2 = sum(hess[(12*n+1):(13*n)])
#                hesssig1 = sum(hess[(18*n+1):(19*n)])
#                hesssig2 = sum(hess[(24*n+1):(25*n)])
               
               # Newton directions
               p = solve(hess)%*%grad
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
          }
          # Calculate std errors
          # Asymptotic variance of theta = inverse of Fisher information matrix
          vals_dl = eval(dl)
          hess = attr(vals_dl,"hessian")
          hess = colSums(hess) # No negative here because negative hessian for negative log-lik-->positive hessian
          theta_var = solve(hess)
          selam = theta_var[1,1]
          semu1 = theta_var[2,2]
          semu2 = theta_var[3,3]
          sesd1 = theta_var[4,4]
          sesd2 = theta_var[5,5]
     }
     
     if (method == "EM") {
          # EM algorithm
          while (diff > tol && iters < maxit) {
               iters = iters + 1
               # E-step
               norm1 = dnorm(y,mu1,sigma1)
               norm1log = dnorm(y,mu1,sigma1,log=TRUE)
               norm2 = dnorm(y,mu2,sigma2)
               norm2log = dnorm(y,mu2,sigma2,log=TRUE)
               condexp_zi = lambda*norm1/(lambda*norm1+(1-lambda)*norm2)
               sumcondexp_zi = sum(condexp_zi)
               Q = sum(condexp_zi*norm1log + (1-condexp_zi)*norm2log)
               # M-step
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
          }
          # Calculate standard errors
          loglikterm = expression(-log((lambda/sigma1)*exp(-(y-mu1)^2/(sigma1^2)) 
                                       + ((1-lambda)/sigma2)*exp(-(y-mu2)^2/(sigma2^2))))
          dl = deriv3(loglikterm,namevec=c("lambda","mu1","mu2","sigma1","sigma2"))
          vals_dl = eval(dl)
          grad = attr(vals_dl,"gradient")
          fisherinfo = t(grad)%*%grad
          theta_var = solve(fisherinfo)
          selam = theta_var[1,1]
          semu1 = theta_var[2,2]
          semu2 = theta_var[3,3]
          sesd1 = theta_var[4,4]
          sesd2 = theta_var[5,5]
     }
     
     return(list(mle=c(lambda=lambda,mu1=mu1,mu2=mu2,sigma1=sigma1,sigma2=sigma2),
                 stderr=c(lambda=selam,mu1=semu1,mu2=semu2,sigma1=sesd1,sigma2=sesd2)))
}