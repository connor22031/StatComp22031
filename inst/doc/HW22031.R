## -----------------------------------------------------------------------------
library(xtable)
library(datasets)
xtable::xtable(head(iris))

## -----------------------------------------------------------------------------
par(mfrow=c(1,1))
plot(lm(Sepal.Length ~ Sepal.Width, data = iris))

## -----------------------------------------------------------------------------
head(iris)

## -----------------------------------------------------------------------------
set.seed(10000)
n=1000
u=runif(n)
a=2
b=2
x=b/((1-u)^(1/a)) #give the derivation

#density histogram of sample
par(mfrow=c(1,1))
hist(x,prob=TRUE,main=bquote(Pareto(2,2)))
y=seq(0,100,0.01)
lines(y,a*(b^a)*1/(y^(a+1)))

## -----------------------------------------------------------------------------
set.seed(10001)
#give the function
samplebeta=function(n,a,b){
  k=0 #counter for accepted
  y=numeric(n)
  while (k<n) {
    u=runif(1)
    x=runif(1) #random variate from g
    if(x^(a-1)*(1-x)^(b-1)>u){
      #we accept x
      k=k+1
      y[k]=x
    }
  }
  y
}

#give the needed sample
x=samplebeta(1000,3,2) 
par(mfrow=c(1,1))
hist(x,prob=TRUE,main=bquote(B(3,2)))
y=seq(0,1,0.01)
lines(y,dbeta(y,3,2))

## -----------------------------------------------------------------------------
set.seed(10002)
#generate a Exp-Gamma mixture
n=1000
r=4
beta=2
lambda=rgamma(n,r,beta) #lambda is random

#now supply the sample of lambda's
x=rexp(n,lambda) #the mixture
x[1:50] #show the first 50 datas to refer to

## -----------------------------------------------------------------------------
#we use the sample from the previous question
par(mfrow=c(1,1))
hist(x,prob=TRUE,main="Gamma-Exp")
y=seq(0,10,0.01)
lines(y,r*(beta^r)*1/((beta+y)^(r+1)))

## -----------------------------------------------------------------------------
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}#one sort
}

## -----------------------------------------------------------------------------
set.seed(12345)
n=c(1e4,2e4,4e4,6e4,8e4)
b=c() #store the time of each simulation
a=c()
for(i in 1:5){
  for(j in 1:100){ #100 simulations
    test=sample(1:n[i])
    b[j]=system.time(quick_sort(test))[1]
  }
  a[i]=mean(b)
}
a

## -----------------------------------------------------------------------------
t=n*log(n)
summary(lm(a~t))#regression
library(ggplot2)
par(mfrow=c(1,1))
ggplot(data=NULL,mapping=aes(x=t,y=a))+geom_point(color="darkred")+geom_smooth(color="blue",formula=y~x,method = "lm")+labs(title="Linear Regression Plot")+theme(plot.title = element_text(hjust = 0.5,size=15))

## -----------------------------------------------------------------------------
set.seed(13245)
m=10000
U1=runif(m)
U2=U1[1:m/2]
T1=exp(U1) #simple MC
T2=(exp(U2)+exp(1-U2))/2 #antithetic variate
mean(T1)
mean(T2)
(var(T1)-var(T2))/var(T1)

## -----------------------------------------------------------------------------
set.seed(10234)
m=10000
theta.hat=se=numeric(2)
g=function(x){
  x^2/sqrt(2*pi)*exp(-x^2/2)*(x>1)
}

u=runif(m) #using f1, inverse transformation method
x=qnorm(u*(1-pnorm(1))+pnorm(1))
fg=g(x)/dnorm(x)*(1-pnorm(1))
theta.hat[1]=mean(fg)
se[1]=sd(fg)

u=runif(m) #using f2, inverse transformation method
x=sqrt(1-2*log(1-u))
fg=g(x)/(x*exp((1-x^2)/2))
theta.hat[2]=mean(fg)
se[2]=sd(fg)

rbind(theta.hat,se) #final estimate and standard deviation
integrate(g,1,Inf)$value #theoretical value

## -----------------------------------------------------------------------------
x=seq(1,5,0.01)
w=2
f1=dnorm(x)*(1-pnorm(1))
f2=x*exp((1-x^2)/2)
g=x^2/sqrt(2*pi)*exp(-x^2/2)

par(mfrow=c(1,1))
plot(x,g,type="l",main="",ylab="",ylim=c(0,1),lwd=w)
lines(x,f1,lty=2,lwd=w)
lines(x,f2,lty=3,lwd=w)
legend("topright",legend=c("g",1:2),lty=1:3,lwd=w,inset=0.02)

## -----------------------------------------------------------------------------
set.seed(10023)
M=10000 #number of replicates
k=5 #number of strata
N=50 #number of time to repeat the estimation
T2=numeric(k)
estimates=numeric(N)
g=function(x){
  exp(-x-log(1+x^2))*(x>0)*(x<1)
}
for(i in 1:N){
  for(j in 1:k){
    u=runif(M/k)
    x=-log(1-(j-1)/k*(1-exp(-1))-u*(1-exp(-1))/k)
    fg=g(x)/(k*exp(-x))*(1-exp(-1))
    T2[j]=mean(fg)
  }
  estimates[i]=mean(T2)
}
mean(estimates)
var(estimates)
k*mean(estimates) #the estimated value

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(123456)

n=20
alpha=0.05
m=1000 #repeated test times
UCL=numeric(m)

for(i in 1:m){
  x=rlnorm(n,0,2)
  y=log(x)
  t1=mean(y)+sd(y)*qt(alpha/2,df=n-1)/sqrt(n)
  t2=mean(y)+sd(y)*qt(1-alpha/2,df=n-1)/sqrt(n)
  if(t1<0 && t2>0) UCL[i]=1 #count the number of intervals that contain mu=0
}

mean(UCL)

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(1234567)

count5test=function(x,y){
  X=x-mean(x)
  Y=y-mean(y)
  outx=sum(X>max(Y))+sum(X<min(Y))
  outy=sum(Y>max(X))+sum(Y<min(X))
  return(as.integer(max(c(outx,outy))>5))
}

m=10000
n1=n2=20
mu1=mu2=0
sigma1=1
sigma2=1.5
alpha=0.055
cfp=fp=numeric(m)

for (i in 1:m) {
  x=rnorm(n1,mu1,sigma1)
  y=rnorm(n2,mu2,sigma2)
  cfp[i]=count5test(x,y)
  fp[i]=var.test(x,y)$p.value
}

mean(cfp) #Count Five test
mean(fp<=0.055) #F test

## -----------------------------------------------------------------------------
set.seed(1234567)
n=seq(10,300,20)
M=length(n)
cfpower=fpower=numeric(M)

for (i in 1:M) {
  n0=n[i]
  for (j in 1:m) {
   x=rnorm(n0,mu1,sigma1)
   y=rnorm(n0,mu2,sigma2)
   cfp[j]=count5test(x,y)
   fp[j]=var.test(x,y)$p.value
  }
 cfpower[i]=mean(cfp) #Count Five test
 fpower[i]=mean(fp<=0.055) #F test
}

par(mfrow=c(1,1))
plot(n,cfpower,ylim=c(0,1),type="l",xlab="n",ylab="power")
lines(n,fpower,lty=2)
legend("bottomright",0.2,c("Count Five","F"),lty=c(1,2),inset=0.02)

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(123456)

t=c(3,5,7,18,43,85,91,98,100,130,230,487)
p=function(x,i){
  1/mean(x[i])
}

library("boot")
boot(t,statistic = p,R=2000)

## -----------------------------------------------------------------------------
set.seed(12344)
m=function(x,i){
  mean(x[i])
}

boot.obj=boot(t,statistic = m,R=2000)
boot.ci(boot.obj,type = c("norm","basic", "perc", "bca"))

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(12343)
n=20
m=1000
me=function(x,i){
  mean(x[i])
}

#calculate the porportion that the intervals miss on the left and right
leftn=rightn=leftb=rightb=leftp=rightp=numeric(m)

for (i in 1:m) {
  x=rnorm(n)
  boot.obj=boot(x,statistic = me,R=2000)
  b=boot.ci(boot.obj,type = c("norm","basic", "perc"))
  if(b$norm[2]>0) leftn[i]=1
  if(b$norm[3]<0) rightn[i]=1
  if(b$basic[4]>0) leftb[i]=1
  if(b$basic[5]<0) rightb[i]=1
  if(b$percent[4]>0) leftp[i]=1
  if(b$percent[5]<0) rightp[i]=1
}

#the empirical coverage rates
1-mean(leftn+rightn)#norm
1-mean(leftb+rightb)#basic
1-mean(leftp+rightp)#percent

#the porportion that the intervals miss

#norm
c(mean(leftn),mean(rightn))
#basic
c(mean(leftb),mean(rightb))
#percent
c(mean(leftp),mean(rightp))

## -----------------------------------------------------------------------------
rm=list(c())
#write the function
pca=function(x,i){
  val=eigen(cov(x[i,]))$values
  return(val[1]/sum(val))
}

library(bootstrap)
n=nrow(scor)
theta.hat <- pca(scor,1:n)
theta.jack <- numeric(n)
for(i in 1:n){
theta.jack[i] <- pca(scor,(1:n)[-i])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
se.jack=se.jack),3)

## -----------------------------------------------------------------------------
rm=list(c())
library(DAAG)
attach(ironslag)
n=length(magnetic)
e1=e2=e3=e4=c()

for (k in 1:(n-1)) {
  for (l in (k+1):n) {
    y=magnetic[-c(k,l)]
    x=chemical[-c(k,l)]
    
    J1=lm(y~x)
    yhat11=J1$coef[1]+J1$coef[2]*chemical[k]
    yhat12=J1$coef[1]+J1$coef[2]*chemical[l]
    e1=c(e1,magnetic[k]-yhat11,magnetic[l]-yhat12)
    
    J2=lm(y~x+I(x^2))
    yhat21=J2$coef[1]+J2$coef[2]*chemical[k]+J2$coef[3]*chemical[k]^2
    yhat22=J2$coef[1]+J2$coef[2]*chemical[l]+J2$coef[3]*chemical[l]^2
    e2=c(e2,magnetic[k]-yhat21,magnetic[l]-yhat22)
    
    J3=lm(log(y)~x)
    yhat31=exp(J3$coef[1]+J3$coef[2]*chemical[k])
    yhat32=exp(J3$coef[1]+J3$coef[2]*chemical[l])
    e3=c(e3,magnetic[k]-yhat31,magnetic[l]-yhat32)
    
    J4=lm(log(y)~log(x))
    yhat41=exp(J4$coef[1]+J4$coef[2]*log(chemical[k]))
    yhat42=exp(J4$coef[1]+J4$coef[2]*log(chemical[l]))
    e4=c(e4,magnetic[k]-yhat41,magnetic[l]-yhat42)
  }
}
c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))
detach(ironslag)

## -----------------------------------------------------------------------------
rm=list(c())
set.seed(11111)
x=rnorm(12)
y=rnorm(12)

R=999
z=c(x,y)
K=1:24
D=numeric(R)
options(warn=-1)
D0=cor(x,y,method = "spearman")
for(i in 1:R){
  k=sample(K,size=12,replace = FALSE)
  x1=z[k]
  y1=z[-k]
  D[i]=cor(x1,y1,method = "spearman")
}
p=mean(c(D0,D)>=D0)
options(warn=0)
p
cor.test(x,y)

## -----------------------------------------------------------------------------
rm=c(list())
set.seed(12321)

#generate the chain
rw.Metropolis=function(sigma,x0,N){
  x=numeric(N)
  x[1]=x0
  u=runif(N)
  k=0
  for(i in 2:N){
    y=rnorm(1,x[i-1],sigma)
    if(u[i]<=exp(abs(x[i-1])-abs(y)))
      x[i]=y else{
        x[i]=x[i-1]
        k=k+1
      }
  }
  return(list(x=x,k=k))
}

N=15000
sigma=c(0.3,0.5,1.5,8)

x0=25
rw1=rw.Metropolis(sigma[1],x0,N)
rw2=rw.Metropolis(sigma[2],x0,N)
rw3=rw.Metropolis(sigma[3],x0,N)
rw4=rw.Metropolis(sigma[4],x0,N)

#the acceptance rates of each chain
print(round(1-c(rw1$k,rw2$k,rw3$k,rw4$k)/N,3))

#compare the chains with different variances
par(mfrow=c(1,1))
rw=cbind(rw1$x,rw2$x,rw3$x,rw4$x)
for(j in 1:4){
  plot(rw[,j],type="l",xlab=bquote(sigma==.(round(sigma[j],3))),ylab="X",ylim=range(rw[,j]))
}
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
}

b=1000 #burn-in length

    #compute diagnostic statistics
    psi <- t(apply(t(rw), 1, cumsum))
    for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))
    print(Gelman.Rubin(psi))

    #plot psi for the four chains
    par(mfrow=c(1,1))
    for (i in 1:4)
        plot(psi[i, (b+1):N], type="l",xlab=i,ylab=bquote(psi))
    par(mfrow=c(1,1)) #restore default
    
    #plot the sequence of R-hat statistics
    rhat <- rep(0, N)
    for (j in (b+1):N)
        rhat[j] <- Gelman.Rubin(psi[,1:j])
    par(mfrow=c(1,1))
    plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
    abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
rm=c(list())
set.seed(123456)

#initialize constants and parameters
N <- 15000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
rho <- .9 #correlation
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2
###### generate the chain #####

X[1, ] <- c(mu1, mu2) #initialize
for (i in 2:N) {
x2 <- X[i-1, 2]
m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
X[i, 1] <- rnorm(1, m1, s1)
x1 <- X[i, 1]
m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
X[i, 2] <- rnorm(1, m2, s2)
}
b <- burn + 1
x <- X[b:N, ]
par(mfrow=c(1,1))
plot(x,main="",cex=.5,xlab=bquote(X[1]),ylab=bquote(X[2]),ylim=range(x[,2]))

L=lm(x[,2]~x[,1])
summary(L)

## -----------------------------------------------------------------------------
par(mfrow=c(1,1))
plot(L$fitted.values,L$residuals)
abline(0,0)
qqnorm(L$residuals)
qqline(L$residuals)
par(mfrow=c(1,1))

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
}

    #compute diagnostic statistics
    psi <- t(apply(t(X), 1, cumsum))
    for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))
    print(Gelman.Rubin(psi))

    #plot psi for the four chains
    par(mfrow=c(1,1))
    for (i in 1:2)
        plot(psi[i, b:N], type="l",xlab=i,ylab=bquote(psi))
    par(mfrow=c(1,1)) #restore default
    
    #plot the sequence of R-hat statistics
    rhat <- rep(0, N)
    for (j in b:N)
        rhat[j] <- Gelman.Rubin(psi[,1:j])
    par(mfrow=c(1,1))
    plot(rhat[b:N], type="l", xlab="", ylab="R")
    abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
rm=c(list())
set.seed(12345)

p=function(N,b1,b2,b3,f0){
  x1=rpois(N,1);x2=rexp(N,1);x3=sample(0:1,N,replace=TRUE)
  g <- function(alpha){
   tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
   p <- 1/(1+tmp)
   mean(p) - f0
  }
  solution <- uniroot(g,c(-20,0))
  return(solution$root)
}

## -----------------------------------------------------------------------------
set.seed(123)
N=1e6;b1=0;b2=1;b3=-1
f0=c(1e-1,1e-2,1e-3,1e-4)
al=c()
for (i in 1:length(f0)) {
  al=c(al,p(N,b1,b2,b3,f0[i]))
}
round(al,3)

## -----------------------------------------------------------------------------
par(mfrow=c(1,1))
plot(log(f0),al)

## -----------------------------------------------------------------------------
rm=c(list())
u=c(11,8,27,13,16,0,23,10,24,2)
v=c(12,9,28,14,17,1,24,11,25,3)

## maximum likelihood
L=function(lambda){
  sum((v-u)/(exp(lambda*(v-u))-1)-u)
}
lambda1=uniroot(L,c(0,0.1))$root
lambda1

## EM
EM=function(u,v,max.it=10000,eps=1e-5){
  lambda=1
  i=1
  lambda1=2
  lambda2=1
  n=length(u)
  x=sum(1/lambda+u-(v-u)/(exp(lambda*(v-u))-1))
  while(abs(lambda1-lambda2)>=eps){
    lambda1=lambda2
    lambda2=n/x
    x=sum(1/lambda2+u-(v-u)/(exp(lambda2*(v-u))-1))
    if(i == max.it) break
    i=i+1
  }
  return(lambda2)
}
EM(u,v,max.it=10000,eps=1e-5)

## -----------------------------------------------------------------------------
a=list(1,2,3,4)
b=unlist(a)
b
is.atomic(b)
c=as.vector(a)
c
is.atomic(c)

## -----------------------------------------------------------------------------
1=='1'
-1<FALSE
'one'<4

## -----------------------------------------------------------------------------
a=c(1,2,3,4)
dim(a)

## -----------------------------------------------------------------------------
x=matrix(c(1,2,3,4,5,6),3,2)
is.matrix(x)
b=array(x)
b
class(b)

## -----------------------------------------------------------------------------
a=data.frame(a=1:2,b=c("a","b"))
as.matrix(a)

## -----------------------------------------------------------------------------
a=matrix(0,nrow=0,ncol=0)
b=as.data.frame(a)
class(b)
nrow(b)
ncol(b)

## -----------------------------------------------------------------------------
scale01=function(x){
  rng=range(x,na.rm = TRUE)
  (x-rng[1])/(rng[2]-rng[1])
}

## -----------------------------------------------------------------------------
lapply(mtcars, scale01)

## -----------------------------------------------------------------------------
lapply(iris[,lapply(iris,class)=="numeric"],scale01)

## -----------------------------------------------------------------------------
vapply(mtcars,sd,0)

## -----------------------------------------------------------------------------
vapply(iris[,vapply(iris,class,'')=='numeric'],sd,0)

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(100)
N=5000
burn=1000
b=burn+1
mu1=mu2=0
sigma1=sigma2=1
rho=0.9

library(Rcpp)
library(StatComp22031)
x1=gibbs_cpp(N,mu1,mu2,sigma1,sigma2,rho)
x2=gibbs_r(N,mu1,mu2,sigma1,sigma2,rho)

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix gibbs_cpp(int N,double mu1,double mu2,double sigma1,double sigma2,double rho){
#    NumericMatrix X(N,2);
#    double s1,s2;
#    s1=sqrt(1-pow(rho,2))*sigma1;
#    s2=sqrt(1-pow(rho,2))*sigma2;
#    X(0,0)=mu1;
#    X(0,1)=mu2;
#    double x1,x2,m1,m2;
#    for(int i=1;i<N;i++){
#      x2=X(i-1,1);
#      m1=mu1+rho*(x2-mu2)*sigma1/sigma2;
#      X(i,0)=rnorm(1,m1,s1)[0];
#      x1=X(i,0);
#      m2=mu2+rho*(x1-mu1)*sigma2/sigma1;
#      X(i,1)=rnorm(1,m2,s2)[0];
#    }
#    return(X);
#  }

## ----eval=FALSE---------------------------------------------------------------
#  gibbs_r=function(N,mu1,mu2,sigma1,sigma2,rho){
#    X=matrix(0,N,2)
#    s1=sqrt(1-rho^2)*sigma1
#    s2=sqrt(1-rho^2)*sigma2
#    X[1,]=c(mu1,mu2)
#    for(i in 2:N){
#      x2=X[i-1,2]
#      m1=mu1+rho*(x2-mu2)*sigma1/sigma2
#      X[i,1]=rnorm(1,m1,s1)
#      x1=X[i,1]
#      m2=mu2+rho*(x1-mu1)*sigma2/sigma1
#      X[i,2]=rnorm(1,m2,s2)
#    }
#    return(X)
#  }

## -----------------------------------------------------------------------------
x1=x1[b:N,]
x2=x2[b:N,]
par(mfrow=c(1,1))
qqplot(x1,x2)

## -----------------------------------------------------------------------------
microbenchmark::microbenchmark(gibbs_cpp(5000,0,0,1,1,0.9),gibbs_r(5000,0,0,1,1,0.9))

