#' @title The value of initial stationary distribution
#' @description Generate a sample of initial stationary distribution
#' @param kappa the parameter of the distribution
#' @param delta the parameter of the distribution
#' @param gamma the parameter of the distribution
#' @param small the truncated variable
#' @return a random sample of size \code{1}
#' @examples
#' \dontrun{
#' sigma0_2=replicate(1000,s0_2(0.5,2,4,0.00001))
#' hist(sigma0_2,freq=FALSE)
#' }
#' @importFrom stats rexp runif
#' @export
s0_2=function(kappa,delta,gamma,small){
  A0=delta*2^kappa*kappa/gamma(1-kappa)
  B=gamma^(1/kappa)/2
  sigma0_2=0
  a0=rexp(1,rate=1)
  decr_u=((a0*kappa)/A0)^(-1/kappa)
  while(decr_u>small){
    e0=rexp(1,rate=B);v0=runif(1)
    sigma0_2=sigma0_2+min(decr_u,e0*(v0)^2)
    a=rexp(1,rate=1);a0=a0+a
    decr_u=((a0*kappa)/A0)^(-1/kappa)
  }
  return(sigma0_2)
}


#' @title Non-compound Poission driving part of TS-OU
#' @description Generate a sample of non-compound Poission driving part
#' @param lambda the parameter of the distribution
#' @param kappa the parameter of the distribution
#' @param delta the parameter of the distribution
#' @param gamma the parameter of the distribution
#' @param Delta the sampling interval
#' @param small the truncated variable
#' @return a random sample of size \code{1}
#' @examples
#' \dontrun{
#' sigma1_2=s1_2(0.1,0.4,2,5,1,0.00001)
#' }
#' @importFrom stats rexp runif
#' @export
s1_2=function(lambda,kappa,delta,gamma,Delta,small){
  A=delta*2^kappa*kappa^2/gamma(1-kappa)
  B=gamma^(1/kappa)/2
  sigma1_2=0
  a1=rexp(1,rate=1)
  decr_u1=((a1*kappa)/(A*lambda*Delta))^(-1/kappa)
  while(decr_u1>small){
    e1=rexp(1,rate=B);v1=runif(1);r1=runif(1)
    sigma1_2=sigma1_2+exp(-lambda*Delta*r1)*min(decr_u1,e1*(v1)^(1/kappa))
    a=rexp(1,rate=1);a1=a1+a
    decr_u1=((a1*kappa)/(A*lambda*Delta))^(-1/kappa)
  }
  return(sigma1_2)
}

#' @title Compound Poission process driving part of TS-OU
#' @description Generate a sample of compound Poission driving part
#' @param lambda the parameter of the distribution
#' @param kappa the parameter of the distribution
#' @param delta the parameter of the distribution
#' @param gamma the parameter of the distribution
#' @param Delta the sampling interval
#' @param small the truncated variable
#' @return a random sample of size \code{1}
#' @examples
#' \dontrun{
#' sigma2_2=s2_2(0.1,0.4,2,5,1,0.00001)
#' }
#' @importFrom stats rexp rgamma
#' @export
s2_2=function(lambda,kappa,delta,gamma,Delta,small){
  t=rexp(1,rate=lambda*delta*gamma*kappa)
  tild_sig=0
  while(t<=Delta){
    re=rgamma(1,1-kappa,rate=1/2*gamma^(1/kappa))
    tild_sig=tild_sig+exp(lambda*t)*re
    t=t+rexp(1,rate=lambda*delta*gamma*kappa)
  }
  sigma2_2=exp(-lambda*Delta)*tild_sig
  return(sigma2_2)
}

#' @title The simulated path of TS-OU process
#' @description Generate a sample of TS-OU process
#' @param lambda the parameter of the distribution
#' @param kappa the parameter of the distribution
#' @param delta the parameter of the distribution
#' @param gamma the parameter of the distribution
#' @param num the number of simulated samples
#' @param Delta the sampling interval
#' @param small the truncated variable
#' @return a random sample of size \code{num}
#' @examples
#' \dontrun{
#' sigma_2=TSOU(0.1,0.4,2,5,1000)
#' plot(sigma_2,type="l")
#' }
#' @export
TSOU=function(lambda,kappa,delta,gamma,num,Delta=1,small=0.00001){
  sigma_2=numeric(num)
  sigma_2[1]=s0_2(kappa,delta,gamma,small)
  for(j in 2:num){
    sigma_2[j]=exp(-lambda*Delta)*sigma_2[j-1]+s1_2(lambda,kappa,delta,gamma,Delta,small)+s2_2(lambda,kappa,delta,gamma,Delta,small)
  }
  return(sigma_2)
}

#' @title A illustration dataset
#' @name data
#' @description A dataset used to illustrate the performance of \code{R_shortest} and \code{Rcpp_shortest}.
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' benchmark(R_shortest(data),Rcpp_shortest(data),
#' columns=c("test","replications","elapsed","relative"))
#' }
NULL

#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{rbenchmark} and \code{microbenchmark} to compare the performance of R functions (\code{R_shortest}) and Cpp functions (\code{Rcpp_shortest}).
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' benchmark(R_shortest(data),Rcpp_shortest(data),
#' columns=c("test","replications","elapsed","relative"))
#' microbenchmark::microbenchmark(gibbs_cpp(5000,0,0,1,1,0.9),
#' gibbs_r(5000,0,0,1,1,0.9))
#' }
#' @import rbenchmark microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm
#' @useDynLib StatComp22031
NULL

#' @title Shortest path model between cities
#' @description Output the connectivity between each city
#' @param M a direct connection situation matrix of size \code{n*2}
#' @return a matrix of size \code{n*n}
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' R_shortest(data)
#' }
#' @export
R_shortest=function(M){
  n=max(M)
  A=matrix(Inf,n,n)
  for(i in 1:nrow(M)){A[M[i,1],M[i,2]]=A[M[i,2],M[i,1]]=1}
  diag(A)=0
  while(TRUE){
    B=A
    for(i in 1:n){
      for(j in 1:n){
        for(k in 1:n){
          if(A[i,j]>A[i,k]+A[k,j]){
            A[i,j]=A[i,k]+A[k,j]
          }
        }
      }
    }
    if(identical(B,A)){break}else{B=A}
  }
  return(A)
}

#' @title Generate texts
#' @name text
#' @description Use R package \code{xtable} to generate examples of texts
#' @examples
#' \dontrun{
#' xtable::xtable(head(iris))
#' }
#' @import xtable
#' @useDynLib StatComp22031
NULL

#' @title Use datasets from R
#' @name dataset
#' @description Use R package \code{datasets} to obtain data
#' @examples
#' \dontrun{
#' head(datasets::iris)
#' }
#' @import datasets
#' @useDynLib StatComp22031
NULL

#' @title Draw pictures for better vision
#' @name graphic
#' @description Use R package \code{graphics} and \code{ggplot2} to draw pictures
#' @examples
#' \dontrun{
#' aes(x = mpg, y = wt)
#' }
#' @import graphics ggplot2
#' @useDynLib StatComp22031
NULL

#' @title The bootstrap method
#' @name bootstrap
#' @description Use R package \code{boot} and \code{bootstrap} to implement bootstrap method
#' @examples
#' \dontrun{
#' t=c(3,5,7,18,43,85,91,98,100,130,230,487)
#' p=function(x,i){
#'  1/mean(x[i])
#' }
#' boot(t,statistic = p,R=2000)
#' }
#' @import boot bootstrap
#' @useDynLib StatComp22031
NULL

#' @title Data Analysis and Graphics Data and Functions
#' @name DAAG
#' @description Use R package \code{DAAG} to get needed data
#' @examples
#' \dontrun{
#' attach(ironslag)
#' }
#' @import DAAG
#' @useDynLib StatComp22031
NULL

#' @title A binary normal Gibbs sampler using R
#' @description Generate a binary normal sample with Gibbs
#' @param N the sample size
#' @param mu1 the mean value of one of the sample
#' @param mu2 the mean value of another sample
#' @param sigma1 the variance of one of the sample
#' @param sigma2 the variance of another sample
#' @param rho the sample correlation coefficient
#' @return A binary normal random sample of size \code{N}
#' @examples
#' \dontrun{
#' gibbs_r(5000,0,0,1,1,0.9)
#' }
#' @export
gibbs_r=function(N,mu1,mu2,sigma1,sigma2,rho){
  X=matrix(0,N,2)
  s1=sqrt(1-rho^2)*sigma1
  s2=sqrt(1-rho^2)*sigma2
  X[1,]=c(mu1,mu2)
  for(i in 2:N){
    x2=X[i-1,2]
    m1=mu1+rho*(x2-mu2)*sigma1/sigma2
    X[i,1]=rnorm(1,m1,s1)
    x1=X[i,1]
    m2=mu2+rho*(x1-mu1)*sigma2/sigma1
    X[i,2]=rnorm(1,m2,s2)
  }
  return(X)
}