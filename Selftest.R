l<-rbinom(1000,1,0.5)
x<-vector(length=1000)

for (i in 1:1000){
  if(l[i]==0){
    x[i]=rnorm(1,mean=1,sd=1)
    
  }
  else{x[i]=rnorm(1,mean=2,sd=2)}
  
}
hist(x)
mixture(x,method="EM",maxit=1000)