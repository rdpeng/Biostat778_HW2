## DER function 
der<-function(xt,y,m){
  gd=matrix(rep(0,5),nrow=1,ncol=5)
  hes=matrix(rep(0,25),nrow=5,ncol=5)
  nor = expression(log(lambda1/sqrt(2*pi)/sqrt(sigma1)*exp((-1)*(y-mu1)^2/(2*sigma1))+(1-lambda1)/sqrt(2*pi)/sqrt(sigma2)*exp((-1)*(y-mu2)^2/(2*sigma2))))
  gra = deriv3(nor,c("lambda1","mu1","mu2","sigma1","sigma2"))
  lambda1=xt[1];mu1=xt[2];mu2=xt[3];sigma1=xt[4];sigma2=xt[5]
  Gra = attr(eval(gra),"gradient")
  Gras = as.matrix(apply(Gra,2,sum))
  ss=matrix(rep(0,25),nrow=5,ncol=5)
  hes = attr(eval(gra),"hessian") 
  Hes = matrix(rep(0,5*5),nrow=5)
  for (i in 1:m){
    #gd=gd+as.matrix(attr(eval(deriv3(~log(lambda1/sqrt(2*pi)/sqrt(sigma1)*exp((-1)*(y-mu1)^2/(2*sigma1))+(1-lambda1)/sqrt(2*pi)/sqrt(sigma2)*exp((-1)*(y-mu2)^2/(2*sigma2))),c("lambda1","mu1","mu2","sigma1","sigma2"))),"gradient"))
    #hes=hes+as.matrix(drop(attr(eval(deriv3(~log(lambda1/sqrt(2*pi)/sqrt(sigma1)*exp((-1)*(y-mu1)^2/(2*sigma1))+(1-lambda1)/sqrt(2*pi)/sqrt(sigma2)*exp((-1)*(y-mu2)^2/(2*sigma2))),c("lambda1","mu1","mu2","sigma1","sigma2"))),"hessian")))
    Hes = Hes + hes[i,,]
    ss=ss+Gra[i,]%*%t(Gra[i,])
  }
  vr=sqrt(diag(solve(ss))); im = -Hes
  IM = sqrt(1/m*diag(solve(im %*% t(im))))
  
  list(gd=Gras,hes=Hes,vr=vr)
}

#
mixture <- function(y, method, maxit = NULL, tol = 1e-08, param0 = NULL) {
  mle=vector(length=5);stderr=vector(length=5)
  if (is.null(param0)==TRUE){
    lambda1=0.5;lambda2=0.5;mu1=10.5;sigma1=59;mu2=20.5;sigma2=249
  }
  else{ lambda1=param0[1];mu1=param0[2];mu2=param0[3];sigma1=param0[4];sigma2=param0[5]
        }
  m<-length(y)
  
  if (method=="EM"){
    it=0
    if (is.null(maxit)==TRUE) {
      maxit=500
    }
    T<-matrix(rep(0,m*2),nrow=m,ncol=2)
    for(i in 1:maxit){
      ll=sum(log(lambda1*dnorm(y,mean=mu1,sd=sqrt(sigma1))+lambda2*dnorm(y,mean=mu2,sd=sqrt(sigma2))))
      ## E
     f1=dnorm(y,mean=mu1,sd=sqrt(sigma1));f2=dnorm(y,mean=mu2,sd=sqrt(sigma2))
      T[,1]=lambda1*f1/((lambda1*f1)+(lambda2*f2))                                            
     T[,2]=lambda2*f2/((lambda1*f1)+(lambda2*f2))   
      ## M
      lambda1=mean(T[,1]);
      lambda2=mean(T[,2])
      mu1=sum(T[1:m,1]*y[1:m])/sum(T[1:m,1]);mu2=sum(T[1:m,2]*y[1:m])/sum(T[1:m,2])
      sigma1=sum(T[1:m,1]*(y[1:m]-mu1)^2)/sum(T[1:m,1]);sigma2=sum(T[1:m,2]*(y[1:m]-mu2)^2)/sum(T[1:m,2])
      ##return
      lln=sum(log(lambda1*dnorm(y,mean=mu1,sd=sqrt(sigma1))+lambda2*dnorm(y,mean=mu2,sd=sqrt(sigma2))))
      if (abs(lln-ll)<= tol){
        break
      }
      it=it+1
     }
     ##STE
   dl=T[,1]/lambda1-T[,2]/lambda2
   dmu1=T[,1]*(y-mu1)/sigma1;dmu2=T[,2]*(y-mu2)/sigma2
   ds1=0.5*T[,1]*((y-mu1)^2-sigma1)/sigma1^2;ds2=0.5*T[,2]*((y-mu2)^2-sigma2)/sigma2^2
   temp=rbind(dl,dmu1,dmu2,ds1,ds2)
  ss=matrix(rep(0,25),nrow=5,ncol=5)
  for (i in 1:m){
    ss=ss+temp[,i]%*%t(temp[,i])
  }
   va=sqrt(diag(solve(ss)))
  mle=c(lambda1,mu1,mu2,sigma1,sigma2);stderr=va
  
  }
  
  ##
  else if (method=="newton"){
    if (is.null(maxit)==TRUE) {
      maxit=100
    }
    x=c(lambda1,mu1,mu2,sigma1,sigma2)
    it=0
    for (i in 1:maxit){
      ll=sum(log(lambda1*dnorm(y,mean=mu1,sd=sqrt(sigma1))+lambda2*dnorm(y,mean=mu2,sd=sqrt(sigma2))))
      gd=der(x,y,m)$gd
      hes=der(x,y,m)$hes
      x=x-solve(hes)%*%gd
      lambda1=x[1];mu1=x[2];mu2=x[3];sigma1=x[4];sigma2=x[5]
      lln=sum(log(lambda1*dnorm(y,mean=mu1,sd=sqrt(sigma1))+lambda2*dnorm(y,mean=mu2,sd=sqrt(sigma2))))
     if (abs(lln-ll)<= tol){
        break
      }
      it=it+1
    }
    lambda1=x[1];mu1=x[2];mu2=x[3];sigma1=x[4];sigma2=x[5]
    va=der(x,y,m)$vr
    mle=c(lambda1,mu1,mu2,sigma1,sigma2);stderr=va
  }
  mle=as.matrix(mle);stderr=as.matrix(stderr)
 rownames(mle)=c("lambda","mu1","mu2","sigma1","sigma2")
  rownames(stderr)=c("lambda","mu1","mu2","sigma1","sigma2")
  ##
  list(mle=mle,stderr=stderr)
}






