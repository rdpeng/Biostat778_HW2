mixture <- function(y, method, maxit = NULL, tol = 1e-08, param0=NULL) {
  Method = match.arg(method,c("Newton","EM"))
  if (is.data.frame(y)){ # transform data.frame to a vector
    y=as.matrix(y)
    y=y[,1]
  }
  
  if (length(param0)==0){ # set initial values using first several stepts of EM
    pk = matrix(c(0.5,-1,5,3,4))
    k=0
    while ( k < 15){
      plast = pk
      loga = log(pk[1])+dnorm(y,mean = pk[2],sd=sqrt(pk[4]),log=T)
      logb = log(1-pk[1])+dnorm(y,mean = pk[3],sd=sqrt(pk[5]),log=T)
      bovera = exp(logb-loga)
      w = (1+bovera)^(-1)
      
      pk[1] = mean(w)
      pk[2] = sum(w*y)/sum(w)
      pk[3] = sum((1-w)*y)/sum(1-w)
      pk[4] = sum(w*(y-plast[2])^2)/sum(w)
      pk[5] = sum((1-w)*(y-plast[3])^2)/sum(1-w)
      k = k + 1
    }
    param0=pk; rm(pk,k)
    # print(param0)
  }
  
  if (Method == "Newton"){
    # define the function that returns function value, gradient and hessian matrix
    func = function(p,data=y,fonly=F){                                                        # p: the parameter vector ; data = y as default
      if (fonly==TRUE)
      {
        if (p[4]>0 & p[5]>0){
          w=p[1];m1=p[2];m2=p[3];v1=p[4];v2=p[5];y=data
          fval = sum(-log(w/sqrt(v1)*dnorm((y-m1)/sqrt(v1))+(1-w)/sqrt(v2)*dnorm((y-m2)/sqrt(v2))))
        }
        else{
          fval=Inf
        }
        result = list(Fval = fval)
      }
      else
      {
        lik = expression(-log(w/sqrt(v1)*dnorm((y-m1)/sqrt(v1))+(1-w)/sqrt(v2)*dnorm((y-m2)/sqrt(v2)))) # negative log likelihood function
        FGH = deriv3(lik,c("w","m1","m2","v1","v2"),function(w,m1,m2,v1,v2,y=data){})
        
        if (p[4]>0 & p[5]>0){
          fgh = FGH(w=p[1],m1=p[2],m2=p[3],v1=p[4],v2=p[5])
          Fval = sum(fgh)
          G = colSums(attr(fgh,"gradient"))
          H = attr(fgh,"hessian")
          k = dim(H)[3]
          Hes = matrix(0,nrow=k,ncol=k)
          for (j in 1:k){
            Hes[j,]=apply(H[,,j],2,sum)
          }
        }
        else {
          Fval = Inf; G=NA; Hes=NA;
        }
        result = list(Fval = Fval, G = matrix(G), H = Hes)
      }
      result
    }
    
    # define a function that modifies the hessian matirx to make it positive definite    
    modNewton = function (H, Beta){
      eig = eigen(H)
      V = eig$vectors
      lambda = eig$values
      n = nrow(H)
      flag = 0
      if (norm(H) == 0){epsilon = 1}
      else {epsilon = norm(H)/Beta}
      
      for (t in 1:n){
        if (lambda[t] <= -epsilon){lambda[t] = -lambda[t]}
        else if (lambda[t] < epsilon){lambda[t] = epsilon}
      }
      
      if (sum(lambda == eig$values)<n){flag = 1}
      B = V %*% diag(lambda) %*% t(V)
      list(B = B, flag = flag)
    }
    
    # define a function that use back tracking to find proper step length    
    armijo = function(fun,p,s,Fval,G){
      t = 0.5; eta = 0.05
      a = 1 ; l = 1
      while (fun(p=p+a[l]*s,fonly=TRUE)$Fval >= Fval + eta*a[l]*t(G)%*%s & l < 30){
        a[l+1] = t * a[l]
        l = l + 1
      }
      alpha = a[length(a)]
      alpha
    }
    
    # define the actual modified Newton function    
    uncMIN = function(fun,p0,maxit,tol){
      pk = as.matrix(p0)
      fgh0 = fun(pk)
      Fval = fgh0$Fval; G = fgh0$G; H = fgh0$H; G0 = fgh0$G
      k = 1 
      while (norm(G)>tol*max(1,norm(G0)) & k<=maxit){
        B = modNewton(H,Beta=10^6)$B
        s = -solve(B,G)                                                  # modified newton direction
        
        alpha = armijo(fun = func, p=pk, s=s, Fval = Fval, G = G)        # the proper step length
        pk_1 = pk + alpha*s                                              
        
        fghk = fun(pk_1)                                                 # compute values for next iter
        Fval = fghk$Fval;G = fghk$G; H = fghk$H
        pk = pk_1
        k=k+1
      }
      param = pk; iter = k-1
      status = as.integer(norm(G)>tol*max(1,norm(G0)))                   # 0-converge, 1-not converge
      info = solve(H)
      se= sqrt(diag(info))
      list(mle = param, stderr = se, Fval = Fval, gradient = G, hessian = H, iter = iter, status = status)
    }
    
    result = uncMIN(fun=func,p0=param0,maxit=maxit,tol=tol)
  }
  else {
    EMmix = function (p0,y,maxit,tol){
      pk = as.matrix(p0)
      plast = pk+10
      k=1
      while (norm(pk-plast) > tol & k < maxit+1){                       
        plast = pk
        loga = log(pk[1])+dnorm(y,mean = pk[2],sd=sqrt(pk[4]),log=T)
        logb = log(1-pk[1])+dnorm(y,mean = pk[3],sd=sqrt(pk[5]),log=T)
        bovera = exp(logb-loga)
        w = (1+bovera)^(-1)     
        # w = pk[1]*dnorm(y,mean = pk[2],sd=sqrt(pk[4]))/(pk[1]*dnorm(y,mean = pk[2],sd=sqrt(pk[4]))+(1-pk[1])*dnorm(y,mean = pk[3],sd=sqrt(pk[5])))
        
        pk[1] = mean(w)
        pk[2] = sum(w*y)/sum(w)
        pk[3] = sum((1-w)*y)/sum(1-w)
        pk[4] = sum(w*(y-plast[2])^2)/sum(w)
        pk[5] = sum((1-w)*(y-plast[3])^2)/sum(1-w)
        k = k + 1
      }
      lik = expression(w*log(lambda/sqrt(v1)*dnorm((y-m1)/sqrt(v1)))+(1-w)*log((1-lambda)/sqrt(v2)*dnorm((y-m2)/sqrt(v2))))
      FGH = deriv3(lik,c("lambda","m1","m2","v1","v2"),function(lambda,m1,m2,v1,v2,w,y){})
      fgh = FGH(lambda=pk[1],m1=pk[2],m2=pk[3],v1=pk[4],v2=pk[5],w=w,y=y)
      Syi = attr(fgh,"gradient")
      SSt = apply(Syi,1,function(x){x%*%t(x)})
      avgSSt = apply(SSt,1,sum)
      Iy = matrix(avgSSt,nrow = 5)
      std.err = sqrt(diag(solve(Iy)))
      iters = k-1
      status = as.integer(norm(pk-plast) > tol)
      list(mle = pk,stderr = std.err, iter = iters, status = status,info = Iy)
    }
    result = EMmix(p0=param0,y=y,maxit = maxit, tol=tol)
  }
  
  result
}
