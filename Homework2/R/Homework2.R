#EM algorithm for mixture of two normal distributions 

mixture <- function(y,method, maxit = NULL, tol = 1e-8,param0=NULL)
  
{ y=as.matrix(y)
  n=length(y)
  
  method=match.arg(method,c("EM","Newton"))
  #####################################################################
  if (method=="EM"){
    ###
    if ( is.null(maxit)){
      if (method=="Newton")
        (maxit=100)
      if (method=="EM")
        (maxit=500)
    }
    if (is.null(param0))
    {th1=c(0.1, 1, 4, 4, 4)}      # theta=(pi1, mu1, mu2, sigma^2, sigma^2), starting value
    else 
    {th1=param0 }
    
    ######
    
    t=vector(length=5)     #t={t[1],t[2],t[3],t[4],t[5])}. where 
    #t[1]=pi1, t[2]=mu1, t[3]=mu2, t[4]=sigma^2, t[5]=sigma^2
    
    y1=function(t) { (t[1]*exp(-((y-t[2])^2)/(2*t[4])))/((2*pi*t[4])^0.5)}   #have 1000 X's
    y2=function(t) { ((1-t[1])*exp(-((y-t[3])^2)/(2*t[5])))/((2*pi*t[5])^0.5)}
    
    z1m=y1(th1)/(y1(th1)+y2(th1))     # z1m=Zij(m)
    z2m=y2(th1)/(y1(th1)+y2(th1))
    
    #sum(z1m)+sum(z2m)    #should add up to n=1000
    
    th2=vector(length=5)
    
    th2[1]=sum(z1m)/(sum(z1m)+sum(z2m))
    th2[2]=(t(z1m)%*%y)/sum(z1m)
    th2[3]=(t(z2m)%*%y)/sum(z2m)
    th2[4]=(t(z1m)%*%((y-th2[2])^2))/sum(z1m)
    th2[5]=(t(z2m)%*%((y-th2[3])^2))/sum(z2m)
    
    th2
    
    count=1
    
    while (max(abs(th2-th1))>tol) {
{th1=th2
 z1m=y1(th1)/(y1(th1)+y2(th1))   
 z2m=y2(th1)/(y1(th1)+y2(th1))
 
 th2=vector(length=5)
 
 th2[1]=sum(z1m)/(sum(z1m)+sum(z2m))
 th2[2]=(t(z1m)%*%y)/sum(z1m)
 th2[3]=(t(z2m)%*%y)/sum(z2m)
 th2[4]=(t(z1m)%*%((y-th2[2])^2))/sum(z1m)
 th2[5]=(t(z2m)%*%((y-th2[3])^2))/sum(z2m)
 if (count>maxit)
   break
}
count=count+1
    }
    
    
    
    
    l <- expression(log(llambda * exp(-(y - lmu1) ^ 2 / (2 * lsigma1)) / sqrt(2 * pi * lsigma1) + 
                          (1 - llambda) * exp(-(y - lmu2) ^ 2 / (2 * lsigma2)) / sqrt(2 * pi * lsigma2)))
    der <- deriv3(l, c('llambda', 'lmu1', 'lmu2', 'lsigma1', 'lsigma2'))
    
    ## mle is the MLE
    llambda <- th2[1]
    lmu1 <- th2[2]
    lmu2 <- th2[3]
    lsigma1 <- th2[4]
    lsigma2 <- th2[5]
    
    grad <- attr(eval(der),"gradient")
    
    score <- matrix(0,5,5)
    for (i in 1:n)
    {
      score <- score + grad[i, ] %*% t(grad[i, ])
    }
    score <- score
    sd <- c(sqrt(diag(solve(score))))
    
    
    
    
    
    x=list(th2[1],th2[2],th2[3],th2[4],th2[5])
    names(x) <- c("lambda","mu1","mu2","sigma1","sigma2")
    print(x)
    
    #list of sd
    y=list(sd[1],sd[2],sd[3],sd[4],sd[5])
    names(y) <- c("lambda","mu1","mu2","sigma1","sigma2")
    print(y)
    
    print(count)
  }                
  
  #END OF EM
  
  
  ############################################################### 
  
  if (method=="Newton"){
    # run EM
    
    ###
    if ( is.null(maxit)){
      if (method=="Newton")
        (maxit=100)
      if (method=="EM")
        (maxit=500)
    }
    if (is.null(param0))
    {th1=c(0.1, 1, 40, 40, 40)}      # theta=(pi1, mu1, mu2, sigma^2, sigma^2), starting value
    else 
    {th1=param0 }
    
    ####
    
    t=vector(length=5)     #t={t[1],t[2],t[3],t[4],t[5])}. where 
    #t[1]=pi1, t[2]=mu1, t[3]=mu2, t[4]=sigma^2, t[5]=sigma^2
    
    y1=function(t) { (t[1]*exp(-((y-t[2])^2)/(2*t[4])))/((2*pi*t[4])^0.5)}   #have 1000 X's
    y2=function(t) { ((1-t[1])*exp(-((y-t[3])^2)/(2*t[5])))/((2*pi*t[5])^0.5)}
    
    z1m=y1(th1)/(y1(th1)+y2(th1))     # z1m=Zij(m)
    z2m=y2(th1)/(y1(th1)+y2(th1))
    
    #sum(z1m)+sum(z2m)    #should add up to n=1000
    
    th2=vector(length=5)
    
    th2[1]=sum(z1m)/(sum(z1m)+sum(z2m))
    th2[2]=(t(z1m)%*%y)/sum(z1m)
    th2[3]=(t(z2m)%*%y)/sum(z2m)
    th2[4]=(t(z1m)%*%((y-th2[2])^2))/sum(z1m)
    th2[5]=(t(z2m)%*%((y-th2[3])^2))/sum(z2m)
    
    th2
    
    count=1
    
    while (max(abs(th2-th1))>tol) {
{th1=th2
 z1m=y1(th1)/(y1(th1)+y2(th1))   
 z2m=y2(th1)/(y1(th1)+y2(th1))
 
 th2=vector(length=5)
 th2[1]=sum(z1m)/(sum(z1m)+sum(z2m))
 th2[2]=(t(z1m)%*%y)/sum(z1m)
 th2[3]=(t(z2m)%*%y)/sum(z2m)
 th2[4]=(t(z1m)%*%((y-th2[2])^2))/sum(z1m)
 th2[5]=(t(z2m)%*%((y-th2[3])^2))/sum(z2m)
 if (count>maxit)
   break
}
count=count+1
    }
    
    #END OF WITHIN EM
    ####################################
    #print(th2)
    
    
    th1[1]=th2[1]
    th1[2]=th2[2]
    th1[3]=th2[3]
    th1[4]=sqrt(th2[4])
    th1[5]=sqrt(th2[5])
    
    th1=th1+0.01
    
    lmix2 <- deriv3(
      ~ -log(p*dnorm((y-u1)/s1)/s1 + (1-p)*dnorm((y-u2)/s2)/s2),
      c("p","u1","u2","s1","s2"),
      function(y,p,u1,u2,s1,s2) NULL)
    
    mix.gr <- function(theta,y){
      p <- theta[1];u1 <- theta[2]; u2 <- theta[3]; s1 <- theta[4]; s2 <- theta[5]
      colSums(attr(lmix2(y,p,u1,u2,s1,s2),"gradient"))}
    
    mix.he <- function(theta,y){
      p <- theta[1];u1 <- theta[2]; u2 <- theta[3]; s1 <- theta[4]; s2 <- theta[5]
      colSums(attr(lmix2(y,p,u1,u2,s1,s2),"hessian"))}
    
    
    th2=th1-solve(mix.he(th1,y))%*%mix.gr(th1,y)
    
    count=1
    
    while (max(abs(th2-th1))>tol) {
{th1=th2
 th2=th1-solve(mix.he(th1,y))%*%mix.gr(th1,y)
 
 if (count>maxit)
   break
}
count=count+1
    }
    
    l <- expression(log(llambda * exp(-(y - lmu1) ^ 2 / (2 * lsigma1)) / sqrt(2 * pi * lsigma1) + 
                          (1 - llambda) * exp(-(y - lmu2) ^ 2 / (2 * lsigma2)) / sqrt(2 * pi * lsigma2)))
    der <- deriv3(l, c('llambda', 'lmu1', 'lmu2', 'lsigma1', 'lsigma2'))
    
    ## mle is the MLE
    llambda <- th2[1]
    lmu1 <- th2[2]
    lmu2 <- th2[3]
    lsigma1 <- th2[4]* th2[4]
    lsigma2 <- th2[5]*th2[5]
    
    grad <- attr(eval(der),"gradient")
    
    score <- matrix(0,5,5)
    for (i in 1:n)
    {
      score <- score + grad[i, ] %*% t(grad[i, ])
    }
    score <- score
    sd <- c(sqrt(diag(solve(score))))
    
    
    
    x=list(th2[1],th2[2],th2[3],th2[4]*th2[4],th2[5]*th2[5])
    names(x) <- c("lambda","mu1","mu2","sigma1","sigma2")
    print(x)
    
    #list of sd
    y=list(sd[1],sd[2],sd[3],sd[4],sd[5])
    names(y) <- c("lambda","mu1","mu2","sigma1","sigma2")
    print(y)
    
    
    
    
    
    
    
    count 
  }  
  ########################################################
  
  # END OF Newton
  
}

