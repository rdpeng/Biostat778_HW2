mixture <- function(y, method, maxit = NULL, tol = 1e-08, param0 = NULL) {
        y=as.vector(unlist(y))
        n=length(y)
        method=match.arg(method,c("EM","newton"))
        
        #Newton's Method
        if(method=="newton"){
                #Set maximum iteration number
                if(!is.numeric(maxit)){
                        maxit=100
                }
                
                #Initialization
                if(is.null(param0)){
                        #To ensure the convergence of Newton's method, I take the advantage of EM method
                        #and set the initialization based on the result of EM method                       
                        emresult=mixture(y,"EM")[[1]]
                        random=rnorm(5,0,0.1)
                        theta=emresult+random
                }else{
                        theta=param0
                }
                
                
                #likelihood function and derivative
                likelihood=expression(log(lambda*exp(-(y-mu1)^2/(2*sigma1))/sqrt(2*pi*sigma1)+(1-lambda)*exp(-(y-mu2)^2/(2*sigma2))/sqrt(2*pi*sigma2)))
                derivative=deriv3(likelihood, c('lambda', 'mu1', 'mu2', 'sigma1', 'sigma2'))
                for(i in 1:maxit){
                        lambda=theta[1]
                        mu1=theta[2]
                        mu2=theta[3]
                        sigma1=theta[4]
                        sigma2=theta[5]
                        #gradient
                        st1=attr(eval(derivative),"gradient")
                        st1=matrix(apply(st1,2,sum),5,byrow=T)
                        #hessian
                        nd2=attr(eval(derivative),"hessian")
                        nd2=data.frame(nd2)
                        nd2=matrix(apply(nd2,2,sum),5,byrow=T)
                        #update parameters
                        theta=solve(nd2,nd2%*%theta-st1)
                        if(sqrt(sum((solve(nd2)%*%st1)^2))<=tol){
                                break;
                        }
                }
                #Information Matrix
                lambda=theta[1]
                mu1=theta[2]
                mu2=theta[3]
                sigma1=theta[4]
                sigma2=theta[5]
                st1=attr(eval(derivative),"gradient")
                score=matrix(0,5,5)
                #########
                for(i in 1:n){
                        score=score+st1[i,]%*%t(st1[i,])
                }
                cov=solve(score)
                std=c(sqrt(diag(cov)))
                return(list(mle=c(lamda=theta[1],mu1=theta[2],mu2=theta[3],sigma1=theta[4],sigma2=theta[5]),stderr=c(lamda=std[1],mu1=std[2],mu2=std[3],sigma1=std[4],sigma2=std[5])))
        }
        
        #EM Method
        if(method=="EM"){    
                #Set maximum iteration number
                if(!is.numeric(maxit)){
                        maxit=500
                }
                #Initialization
                e=numeric(0)
                lambdaem=numeric(0)
                mu1em=numeric(0)
                mu2em=numeric(0)
                sigma1em=numeric(0)
                sigma2em=numeric(0)
                if(is.null(param0)){
                        lambdaem[1]=0.1
                        mu1em[1]=2
                        mu2em[1]=3
                        sigma1em[1]=4
                        sigma2em[1]=5
                }else{
                        lambdaem[1]=param0[1]
                        mu1em[1]=param0[2]
                        mu2em[1]=param0[3]
                        sigma1em[1]=param0[4]
                        sigma2em[1]=param0[5]
                }

                thetaem=matrix(0,maxit,5)
                for(j in 1:maxit){                       
                        e=lambdaem[j]*dnorm(y,mu1em[j],sigma1em[j]^(1/2))/(lambdaem[j]*dnorm(y,mu1em[j],sigma1em[j]^(1/2))+(1-lambdaem[j])*dnorm(y,mu2em[j],sigma2em[j]^(1/2)))        
                        lambdaem[j+1]=mean(e)
                        mu1em[j+1]=sum(y*e)/sum(e)
                        mu2em[j+1]=sum(y*(1-e))/sum(1-e)
                        sigma1em[j+1]=sum(e*((y-mu1em[j+1])^2))/sum(e)
                        sigma2em[j+1]=sum((1-e)*((y-mu2em[j+1])^2))/sum(1-e)
                        #tolerance check,use norm
                        if(sqrt((lambdaem[j+1]-lambdaem[j])^2+(mu1em[j+1]-mu1em[j])^2+(mu2em[j+1]-mu2em[j])^2+(sigma1em[j+1]-sigma1em[j])^2+(sigma2em[j+1]-sigma2em[j])^2)<=tol){
                                break;
                        }        
                }
                lambda=lambdaem[j+1]
                mu1=mu1em[j+1]
                mu2=mu2em[j+1]
                sigma1=sigma1em[j+1]
                sigma2=sigma2em[j+1]
                
                #Information Matrix
                likelihood=expression(log(lambda*exp(-(y-mu1)^2/(2*sigma1))/sqrt(2*pi*sigma1)+(1-lambda)*exp(-(y-mu2)^2/(2*sigma2))/sqrt(2*pi*sigma2)))
                derivative=deriv3(likelihood, c('lambda', 'mu1', 'mu2', 'sigma1', 'sigma2'))

                lambda=lambda
                mu1=mu1
                mu2=mu2
                sigma1=sigma1
                sigma2=sigma2
                st1=attr(eval(derivative),"gradient")
                score=matrix(0,5,5)
                for(i in 1:n){
                        score=score+st1[i,]%*%t(st1[i,])
                }

                cov=solve(score)
                se=c(sqrt(diag(cov)))
                return(list(mle=c(lambda=lambda,mu1=mu1,mu2=mu2,sigma1=sigma1,sigma2=sigma2),stderr=c(lamda=se[1],mu1=se[2],mu2=se[3],sigma1=se[4],sigma2=se[5])))
        }        
}