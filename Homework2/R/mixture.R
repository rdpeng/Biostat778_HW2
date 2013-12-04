mixture <- function(y, method, maxit = NULL, tol = 1e-08, param0 = NULL) {
        y=unlist(y)
        n=length(y)
        method=match.arg(method,c("EM","newton"))
        
        ##EM Algorithm
        if(method=="EM"){
                
                ##Starting value
                e=numeric(0)
                tlam=numeric(0)
                tmu1=numeric(0)
                tmu2=numeric(0)
                tsigma1=numeric(0)
                tsigma2=numeric(0)
                if(is.null(param0)){
                        tlam[1]=0.1
                        tmu1[1]=3
                        tmu2[1]=2
                        tsigma1[1]=4
                        tsigma2[1]=5
                }else{
                        tlam[1]=param0[1]
                        tmu1[1]=param0[2]
                        tmu2[1]=param0[3]
                        tsigma1[1]=param0[4]
                        tsigma2[1]=param0[5]
                }
                if(!is.numeric(maxit)){
                        maxit=500
                }
                for(j in 1:maxit){
                        
                        e=tlam[j]*dnorm(y,tmu1[j],tsigma1[j]^(1/2))/(tlam[j]*dnorm(y,tmu1[j],tsigma1[j]^(1/2))+(1-tlam[j])*dnorm(y,tmu2[j],tsigma2[j]^(1/2)))        
                        tlam[j+1]=mean(e)
                        tmu1[j+1]=sum(y*e)/sum(e)
                        tmu2[j+1]=sum(y*(1-e))/sum(1-e)
                        tsigma1[j+1]=sum(e*((y-tmu1[j+1])^2))/sum(e)
                        tsigma2[j+1]=sum((1-e)*((y-tmu2[j+1])^2))/sum(1-e)
                        if(sqrt((tlam[j+1]-tlam[j])^2+(tmu1[j+1]-tmu1[j])^2+(tmu2[j+1]-tmu2[j])^2+(tsigma1[j+1]-tsigma1[j])^2+(tsigma2[j+1]-tsigma2[j])^2)<=tol){
                                break;
                        }        
                }
                lambda=tlam[j+1]
                mu1=tmu1[j+1]
                mu2=tmu2[j+1]
                sigma1=tsigma1[j+1]
                sigma2=tsigma2[j+1]
                
                ##Information Matrix
                #d=deriv3(~log(lambda*exp(-(y-mu1)^2/(2*sigma1))/sqrt(2*pi*sigma1)+(1-lambda)*exp(-(y-mu2)^2/(2*sigma2))/sqrt(2*pi*sigma2)),c("lambda","mu1","mu2","sigma1","sigma2"),function.arg=TRUE)
                l=expression(log(lambda * exp(-(y - mu1) ^ 2 / (2 * sigma1)) / sqrt(2 * pi * sigma1) + (1 - lambda) * exp(-(y - mu2) ^ 2 / (2 * sigma2)) / sqrt(2 * pi * sigma2)))
                de=deriv3(l, c('lambda', 'mu1', 'mu2', 'sigma1', 'sigma2'))
                
                lambda=lambda
                mu1=mu1
                mu2=mu2
                sigma1=sigma1
                sigma2=sigma2
                grad=attr(eval(de),"gradient")
                score=matrix(0,5,5)
                for(i in 1:n){
                        score=score+grad[i,]%*%t(grad[i,])
                }
                #score=score*(1/n)
                cov=solve(score)
                se=c(sqrt(diag(cov)))
                return(list(mle=c(lambda=lambda,mu1=mu1,mu2=mu2,sigma1=sigma1,sigma2=sigma2),stderr=c(lamda=se[1],mu1=se[2],mu2=se[3],sigma1=se[4],sigma2=se[5])))
        }
        
        
        ##Newton Method
        if(method=="newton"){
                
                if(is.null(param0)){
                        temp=mixture(y,"EM")[[1]]
                        add=runif(4,-0.5,0.5)
                        add=as.vector(cbind(runif(1,-0.05,0.05),t(add)))
                        theta=temp+add
                        #theta=c(0.5,10.5,20.8,60,250)
                }else{
                        theta=param0
                }
                if(!is.numeric(maxit)){
                        maxit=100
                }
                l=expression(log(lambda * exp(-(y - mu1) ^ 2 / (2 * sigma1)) / sqrt(2 * pi * sigma1) + (1 - lambda) * exp(-(y - mu2) ^ 2 / (2 * sigma2)) / sqrt(2 * pi * sigma2)))
                #d=deriv3(~log(lambda*exp(-(y-mu1)^2/(2*sigma1))/sqrt(2*pi*sigma1)+(1-lambda)*exp(-(y-mu2)^2/(2*sigma2))/sqrt(2*pi*sigma2)),c("lambda","mu1","mu2","sigma1","sigma2"),function.arg=TRUE)
                de=deriv3(l, c('lambda', 'mu1', 'mu2', 'sigma1', 'sigma2'))
                for(i in 1:maxit){
                        lambda=theta[1]
                        mu1=theta[2]
                        mu2=theta[3]
                        sigma1=theta[4]
                        sigma2=theta[5]
                        
                        grad=attr(eval(de),"gradient")
                        grad=matrix(apply(grad,2,sum),5,byrow=T)
                        hes=attr(eval(de),"hessian")
                        hes=data.frame(hes)
                        hes=matrix(apply(hes,2,sum),5,byrow=T)
                        theta=theta-solve(hes)%*%grad
                        if(sqrt(sum((solve(hes)%*%grad)^2))<=tol){
                                break;
                        }
                }
                ##Information Matrix
                lambda=theta[1]
                mu1=theta[2]
                mu2=theta[3]
                sigma1=theta[4]
                sigma2=theta[5]
                grad=attr(eval(de),"gradient")
                score=matrix(0,5,5)
                for(i in 1:n){
                        score=score+grad[i,]%*%t(grad[i,])
                }
                #score=score*(1/n)
                cov=solve(score)
                se=c(sqrt(diag(cov)))
                return(list(mle=c(lamda=theta[1],mu1=theta[2],mu2=theta[3],sigma1=theta[4],sigma2=theta[5]),stderr=c(lamda=se[1],mu1=se[2],mu2=se[3],sigma1=se[4],sigma2=se[5])))
        }
        
}