mixture <-
function(y, method=c("newton","EM"), maxit = NULL, tol = 1e-08, param0 = NULL)
{
y = unlist(y)
n = length(y)
method=match.arg(method,c("EM","newton"))
# Newton method
if (method == "newton")
{
if (is.null(maxit)) maxit = 100
if (is.null(param0))
theta=c(0.5,10.5,20.8,60,250)
else
theta=param0
beta = expression(log(lambda/sqrt(2*pi*sig1)*exp(-(y-mu1)^2/(2*sig1))+(1-lambda)/sqrt(2*pi*sig2)*exp(-(y-mu2)^2/(2*sig2))))
gr = deriv3(beta,c("lambda","mu1","mu2","sig1","sig2"))
for (i in 1:maxit)
{
lambda=theta[1]
mu1=theta[2]
sig1=theta[4]
mu2=theta[3]
sig2=theta[5]
grad = attr(eval(gr),"gradient")
Grad = as.matrix(apply(grad,2,sum))
hes = attr(eval(gr),"hessian")
Hes = matrix(rep(0,5*5),nrow=5)
for(j in 1:n)
{
Hes = Hes + hes[j,,]
}
theta = theta - solve(Hes) %*% Grad
# tolerance
tolr = (lambda-theta[1])^2+(mu1-theta[2])^2+(mu2-theta[3])^2+(sig1-theta[4])^2+(sig2-theta[5])^2
if(tolr < 1e-08)
break
}
lambda=theta[1]
mu1=theta[2]
sig1=theta[4]
mu2=theta[3]
sig2=theta[5]
# Information Matrix
grad = attr(eval(gr),"gradient")
SS = matrix(rep(0,25),ncol=5)
for(i in 1:n){
SS=SS+grad[i,]%*%t(grad[i,])
}
SE = c(sqrt(diag(solve(SS))))
}
else if (method == "EM")
{
if (is.null(maxit))  maxit = 500
if (is.null(param0))
{
lambda = 0.5
mu1 = 100
mu2 = 2
sig1 = 60
sig2 = 150
}
else
{
lambda = param0[1]
mu1 = param0[2]
mu2 = param0[3]
sig1 = param0[4]
sig2 = param0[5]
}
Tmax = matrix(rep(0,2*n),ncol=2)
for(i in 1:maxit)
{
ilambda=lambda
imu1=mu1
isig1=sig1
imu2=mu2
isig2=sig2
# E step
f1=dnorm(y,mu1,sqrt(sig1))
f2=dnorm(y,mu2,sqrt(sig2))
Tmax[,1] = lambda*f1/(lambda*f1+(1-lambda)*f2)
Tmax[,2] = (1-lambda)*f2/(lambda*f1+(1-lambda)*f2)
# M step
lambda = sum(Tmax[,1])/n
mu1 = sum(Tmax[,1]*y)/sum(Tmax[,1])
mu2 = sum(Tmax[,2]*y)/sum(Tmax[,2])
sig1 = sum(Tmax[,1]*(y-mu1)^2)/sum(Tmax[,1])
sig2 = sum(Tmax[,2]*(y-mu2)^2)/sum(Tmax[,2])
# tolerance
tolr = (lambda-ilambda)^2+(mu1-imu1)^2+(mu2-imu2)^2+(sig1-isig1)^2+(sig2-isig2)^2
if(tolr < 1e-08)
break
}
# Information Matrix
bp = lambda*dnorm(y,mu1,sqrt(sig1))/(lambda*dnorm(y,mu1,sqrt(sig1))+(1-lambda)*dnorm(y,mu2,sqrt(sig2)))
Slambda = bp/lambda-(1-bp)/(1-lambda)
Smu1 = bp*(y-mu1)/sig1^2
Smu2 = (1-bp)*(y-mu2)/sig2^2
Ssig1 = bp/(2*sig1)*((y-mu1)^2/sig1-1)
Ssig2 = (1-bp)/(2*sig2)*((y-mu2)^2/sig2-1)
S = rbind(Slambda,Smu1,Smu2,Ssig1,Ssig2)
SS = matrix(rep(0,25),ncol=5)
for(i in 1:n)
{
Si = (S[,i])
SS = SS + Si %*% t(Si)
}
SE = c(sqrt(diag(solve(SS))))
}
# Output
list(mle=c(lambda=lambda,mu1=mu1,mu2=mu2,sigma1=sig1,sigma2=sig2),stderr=c(lambda=SE[1],mu1=SE[2],mu2=SE[3],sigma1=SE[4],sigma2=SE[5]))
}
