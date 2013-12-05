#' Mixture models function
#'
#' \code{mixture} Fits a mixture model with either EM or newton
#'
#' @param y data vector
#' @param method either 'EM' or 'newton'
#' @param maxit number of iterations for either algorithm
#' @param tol how much any parameter can differ from the previous fitted parameter. If any parameter differs by more than this much, we keep searching
#' @param param0 A starting value for the MLE. If none are provided, we calculate a nieve estimate using the quantiles of y.
mixture <- function(y, method='EM', maxit = 10000, tol = 1e-08, param0 = NULL) {

match.arg(method,choices=c('EM','newton'))
if(is.null(maxit)) maxit<-1000
#initialize:
#I converted everything to variance, instead of std.dev.
if(!is.null(param0)&'sigma1'%in%names(param0)){
	param0[['sigma1']]<-param0[['sigma1']]^2
	param0[['sigma2']]<-param0[['sigma2']]^2
	names(param0)<-c('lambda','mu1','mu2','var1','var2')
}#Above code is sloppy, but it works
n <- length(y)
y <- y[order(y)] #for creating initial groups


#If param0 is null
lambda0<-.5 #for use if param0 is null
Z0<-rep(FALSE,n)
Z0[1:floor(n*lambda0)]<-TRUE
if(is.null(param0)){
	param0<-list(lambda=lambda0,mu1=mean(y[Z0]),mu2=mean(y[!Z0]),var1=var(y[Z0]),var2=var(y[!Z0]))
}
param<-param0



#################################
#################################
#################################
#EM METHOD
if(method=='EM'){

iter<-0
tolNotMet<-TRUE#tolerance stop
while(iter<maxit&tolNotMet){
	lastParam<-param

	#updated E(Z)
	EZ_numerator <-  param$lambda*dnorm(y,mean=param$mu1,sd=sqrt(param$var1))
	EZ_denominator<-EZ_numerator+(1-param$lambda)*dnorm(y,mean=param$mu2,sd=sqrt(param$var2))
	EZ<- EZ_numerator/EZ_denominator

	#Update parameter estimates
	param$lambda<-mean(EZ)
	sumEZ<-sum(EZ)
	param$mu1<-sum(EZ*y)/sumEZ
	param$mu2<-sum((1-EZ)*y)/(n-sumEZ)
	param$var1<-sum(EZ*(y-param$mu1)^2)/sumEZ
	param$var2<-sum((1-EZ)*(y-param$mu2)^2)/(n-sumEZ)

	iter<-iter+1
	if(abs(max(unlist(lastParam)-unlist(param)))<tol) tolNotMet<-FALSE
} #could then recalculate EZ, but the tolerance parameter should ensure that it doesn't change anyway

#Get expected full data score, conditional on observed data
Escore_yizi<-matrix(NA,n,5)
colnames(Escore_yizi)<-c('mu1','mu2','lambda','sigma1','sigma2')

Escore_yizi[,'mu1']<-EZ*(y-param$mu1)/(param$var1)
Escore_yizi[,'mu2']<-(1-EZ)*(y-param$mu2)/(param$var2)
Escore_yizi[,'sigma1']<- EZ*(-1/sqrt(param$var1) + (y-param$mu1)^2/sqrt(param$var1)^3)
Escore_yizi[,'sigma2']<- (1-EZ)*(-1/sqrt(param$var2) + (y-param$mu2)^2/sqrt(param$var2)^3)
Escore_yizi[,'lambda']<- EZ/param$lambda -(1-EZ)/(1-param$lambda)

sum_Escore_yizi<-crossprod(Escore_yizi)
mle_stderr<-sqrt(diag(solve(sum_Escore_yizi)))


}################################
#################################
#################################






#################################
#################################
#################################
#NEWTON's METHOD
if(method=='newton'){
#define nlly, negative log likelihood of y
#we'll take the derivatives of this, and then sum the derivatives
nlly<- 
~ -log(
lambda/sqrt(2*pi*var1) * exp( ((y-mu1)^2) / (-2*var1) )  +
(1-lambda)/sqrt(2*pi*var2) * exp( ((y-mu2)^2) / (-2*var2) )
) 
#end def
deriv3_over_y<-deriv3(nlly,names(param0)) 


iter<-0
tolNotMet<-TRUE#tolerance stop
while(iter<maxit&tolNotMet){#NEED TO ADD TOLERENCE PART!!!!
	# Re-assign variables before evaluating derivatives for each y_i
	lastParam<-param

	for(name in names(param)) assign(name,param[[name]])
	evalderiv3<-eval(deriv3_over_y)
	gradient<-colSums(attr(evalderiv3,'gradient'))
	hessian<-apply(attr(evalderiv3,'hessian'),2:3,sum)

	nextDirection<- solve(hessian)%*%gradient
	param<-as.list(unlist(param)+nextDirection)
	names(param)<-names(param0)

	#error messages
	if(param$lambda<0|param$lambda>1) stop(paste0(iter,'th iteration gives a lambda value outside [0,1]'))
	if(param$var1<0)stop(paste0(iter,'th iteration gives a negative variance for group 1'))
	if(param$var2<0)stop(paste0(iter,'th iteration gives a negative variance for group 2'))

	iter<-iter+1
	if(abs(max(unlist(lastParam)-unlist(param)))<tol) tolNotMet<-FALSE
}#End Loop


solveInfo<-solve(-hessian)
#last warning
if(any(diag(solveInfo)<0)) warning('Inverse of covariance matrix has negative diagonal elements, which is not allowed for variance values.')
mle_stderr<-sqrt(diag(solveInfo))
}

#################################
#################################
#################################
#Get output


mle<-list()
mle[['lambda']]<-param[['lambda']]
mle[['mu1']]<-param[['mu1']]
mle[['mu2']]<-param[['mu2']]
mle[['sigma1']]<-sqrt(param[['var1']])
mle[['sigma2']]<-sqrt(param[['var2']])

return(list(mle=mle,stderr=mle_stderr,iter=iter))

## Return a list with elements `mle' for the maximum likelhood estimates and
## `stderr' for their standard errors.
}

