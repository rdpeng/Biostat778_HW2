der <-
function(xt,y,m){
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
