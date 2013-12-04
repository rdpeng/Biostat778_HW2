\name{mixture}
\alias{mixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute parameters for mixture of 2 normal distributions
}
\description{
Takes in a vector of values and fits a mixture of 2 normal distributions to the data. Outputs MLEs and their standard errors, using Newton or EM algorithm. 
}
\usage{
mixture(y, method, maxit = NULL, tol = 1e-08, param0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A vector of values to fit the mixture distribution to
}
  \item{method}{
One of either "Newton" or "EM"
}
  \item{maxit}{
Maximum number of iterations. If NULL, 100 for Newton method and 500 for EM method
}
  \item{tol}{
Tolerance for convergence. Default value is 1E-8
}
  \item{param0}{
Vector of initial estimates in the following form: c(lambda,mu1,mu2,sigma1,sigma2) where sigma1 and sigma2 are the variances, not the standard deviations. 
}
}
\details{
None.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
	\item{mle }{Vector of MLEs for lambda, mu1, mu2, sigma1, and sigma2}
	\item{stderr }{Vector of standard errors for MLEs of lambda, mu1, mu2, sigma1, and sigma2}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Shan Andrews
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
y1<-rnorm(1000,10,5)
y2<-rnorm(1000,20,10)

## The function is currently defined as
mixture <- function(y, method, maxit = NULL, tol = 1e-08, param0 = NULL) {

method<-match.arg(method,c("Newton","newton","NEWTON","EM","em")) 
method<-ifelse(nchar(method)>2,"Newton","EM")  
maxit<-ifelse(is.null(maxit),ifelse(method=="Newton",100,500),maxit)
if (method=="Newton"){
#Not newton, but put EM here so it wouldn't error
#collect<-matrix(0,2,5)
#for (i in 1:100){
##Come up with ways to initialize mu and sigma values: cluster y into 2 groups, seed with extrema
if(length(param0)>0){
lambda<-param0[1]
mu1.i<-param0[2]
mu2.i<-param0[3]
sigma1<-param0[4]
sigma2<-param0[5]
}
if(length(param0)==0){
	fit<-kmeans(y,centers=c(min(y),max(y)))
	clust1<-y[which(fit$cluster==1)]
	clust2<-y[which(fit$cluster==2)]
mu1.i<-mean(clust1)
sigma1<-var(clust1)
mu2.i<-mean(clust2)
sigma2<-var(clust2)
lambda<-length(which(fit$cluster==1))/length(y)}

	
count<-1
gap<-1000

dist1<-(lambda/sqrt(sigma1))*exp((-(y-mu1.i)^2)/(2*sigma1))
dist2<-((1-lambda)/sqrt(sigma2))*exp((-(y-mu2.i)^2)/(2*sigma2))

temp<-sum(log(dist1+dist2))
while((gap>tol)&(count<maxit)){
lambda<-mean(dist1/(dist1+dist2))
mu1<-mean((dist1/(dist1+dist2))*y)/lambda	
mu2<-(mean(y)-lambda*mu1)/(1-lambda)

sigma1<-mean((dist1/(dist1+dist2))*((y-mu1.i)^2))/lambda
sigma2<-mean((1-(dist1/(dist1+dist2)))*(y-mu2.i)^2)/(1-lambda)
mu1.i<-mu1
mu2.i<-mu2

dist1<-(lambda/sqrt(sigma1))*exp((-(y-mu1)^2)/(2*sigma1))
dist2<-((1-lambda)/sqrt(sigma2))*exp((-(y-mu2)^2)/(2*sigma2))

here<-temp
temp<-sum(log(dist1+dist2))
gap<-abs(temp-here)
count<-count+1
}
#collect[i,]<-c(lambda,mu1,mu2,sigma1,sigma2)
#}
# lambda.se<-sqrt(var(collect[,1])/100)
# mu1.se<-sqrt(var(collect[,2])/100)
# mu2.se<-sqrt(var(collect[,3])/100)
# sigma1.se<-sqrt(var(collect[,4])/100)
# sigma2.se<-sqrt(var(collect[,5])/100)
lambda.se<-((1/lambda)+(1/(1-lambda)))/length(y)
mu1.se<-abs(mean((y-mu1)/sigma1))
mu2.se<-abs(mean((y-mu2)/sigma2))
sigma1.se<-abs(mean((-1/sqrt(sigma1))+(1/sqrt(sigma1)^3)*(y-mu1)))
sigma2.se<-abs(mean((-1/sqrt(sigma2))+(1/sqrt(sigma2)^3)*(y-mu2)))
}
if (method=="EM"){
#collect<-matrix(0,2,5)
#for (i in 1:100){
##Come up with ways to initialize mu and sigma values: cluster y into 2 groups, seed with extrema
if(length(param0)>0){
lambda<-param0[1]
mu1.i<-param0[2]
mu2.i<-param0[3]
sigma1<-param0[4]
sigma2<-param0[5]
}
if(length(param0)==0){
	fit<-kmeans(y,centers=c(min(y),max(y)))
	clust1<-y[which(fit$cluster==1)]
	clust2<-y[which(fit$cluster==2)]
mu1.i<-mean(clust1)
sigma1<-var(clust1)
mu2.i<-mean(clust2)
sigma2<-var(clust2)
lambda<-length(which(fit$cluster==1))/length(y)}

	
count<-1
gap<-1000

dist1<-(lambda/sqrt(sigma1))*exp((-(y-mu1.i)^2)/(2*sigma1))
dist2<-((1-lambda)/sqrt(sigma2))*exp((-(y-mu2.i)^2)/(2*sigma2))

temp<-sum(log(dist1+dist2))
while((gap>tol)&(count<maxit)){
lambda<-mean(dist1/(dist1+dist2))
mu1<-mean((dist1/(dist1+dist2))*y)/lambda	
mu2<-(mean(y)-lambda*mu1)/(1-lambda)

sigma1<-mean((dist1/(dist1+dist2))*((y-mu1.i)^2))/lambda
sigma2<-mean((1-(dist1/(dist1+dist2)))*(y-mu2.i)^2)/(1-lambda)
mu1.i<-mu1
mu2.i<-mu2

dist1<-(lambda/sqrt(sigma1))*exp((-(y-mu1)^2)/(2*sigma1))
dist2<-((1-lambda)/sqrt(sigma2))*exp((-(y-mu2)^2)/(2*sigma2))

here<-temp
temp<-sum(log(dist1+dist2))
gap<-abs(temp-here)
count<-count+1
}
#collect[i,]<-c(lambda,mu1,mu2,sigma1,sigma2)
#}
# lambda.se<-sqrt(var(collect[,1])/100)
# mu1.se<-sqrt(var(collect[,2])/100)
# mu2.se<-sqrt(var(collect[,3])/100)
# sigma1.se<-sqrt(var(collect[,4])/100)
# sigma2.se<-sqrt(var(collect[,5])/100)
lambda.se<-((1/lambda)+(1/(1-lambda)))/length(y)
mu1.se<-abs(mean((y-mu1)/sigma1))
mu2.se<-abs(mean((y-mu2)/sigma2))
sigma1.se<-abs(mean((-1/sqrt(sigma1))+(1/sqrt(sigma1)^3)*(y-mu1)))
sigma2.se<-abs(mean((-1/sqrt(sigma2))+(1/sqrt(sigma2)^3)*(y-mu2)))
}
output<-c(lambda,mu1,mu2,sqrt(sigma1),sqrt(sigma2))
names(output)<-c("lambda","mu1","mu2","sigma1","sigma2")
ses<-c(lambda.se,mu1.se,mu2.se,sigma1.se,sigma2.se)
names(ses)<-c("lambda","mu1","mu2","sigma1","sigma2")
mylist<-list("mle"=output,"stderr"=ses)
return(mylist)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line