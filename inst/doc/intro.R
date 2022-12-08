## ----eval=FALSE---------------------------------------------------------------
#  function(kappa,delta,gamma,small){
#    A0=delta*2^kappa*kappa/gamma(1-kappa)
#    B=gamma^(1/kappa)/2
#    sigma0_2=0
#    a0=rexp(1,rate=1)
#    decr_u=((a0*kappa)/A0)^(-1/kappa)
#    while(decr_u>small){  ##control truncation number
#      e0=rexp(1,rate=B);v0=runif(1)
#      sigma0_2=sigma0_2+min(decr_u,e0*(v0)^2)
#      a=rexp(1,rate=1);a0=a0+a
#      decr_u=((a0*kappa)/A0)^(-1/kappa)
#    }
#    return(sigma0_2)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp22031)
set.seed(50)
sigma0_2=replicate(1000,s0_2(0.5,2,4,0.00001))
par(mfrow=c(1,1))
hist(sigma0_2,freq=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  function(lambda,kappa,delta,gamma,Delta,small){
#    A=delta*2^kappa*kappa^2/gamma(1-kappa)
#    B=gamma^(1/kappa)/2
#    sigma1_2=0
#    a1=rexp(1,rate=1)
#    decr_u1=((a1*kappa)/(A*lambda*Delta))^(-1/kappa)
#    while(decr_u1>small){
#      e1=rexp(1,rate=B);v1=runif(1);r1=runif(1)
#      sigma1_2=sigma1_2+exp(-lambda*Delta*r1)*min(decr_u1,e1*(v1)^(1/kappa))
#      a=rexp(1,rate=1);a1=a1+a
#      decr_u1=((a1*kappa)/(A*lambda*Delta))^(-1/kappa)
#    }
#    return(sigma1_2)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  function(lambda,kappa,delta,gamma,Delta,small){
#    t=rexp(1,rate=lambda*delta*gamma*kappa)
#    tild_sig=0
#    while(t<=Delta){
#      re=rgamma(1,1-kappa,rate=1/2*gamma^(1/kappa))
#      tild_sig=tild_sig+exp(lambda*t)*re
#      t=t+rexp(1,rate=lambda*delta*gamma*kappa)
#    }
#    sigma2_2=exp(-lambda*Delta)*tild_sig
#    return(sigma2_2)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  function(lambda,kappa,delta,gamma,num,Delta=1,small=0.00001){
#    sigma_2=numeric(num)
#    sigma_2[1]=s0_2(kappa,delta,gamma,small)
#    for(j in 2:num){
#      sigma_2[j]=exp(-lambda*Delta)*sigma_2[j-1]+s1_2(lambda,kappa,delta,gamma,Delta,small)+s2_2(lambda,kappa,delta,gamma,Delta,small)
#    }
#    return(sigma_2)
#  }

## ----eval=TRUE----------------------------------------------------------------
set.seed(10)
sigma_2=TSOU(0.02,0.4,2,5,5000)
r=function(lambda,s){exp(-lambda*abs(s))}
x=0:35
par(mfrow=c(1,1))
acf(sigma_2)
lines(x,r(0.02,x),type="l")

## ----eval=FALSE---------------------------------------------------------------
#  function(M){
#    n=max(M)
#    A=matrix(Inf,n,n)
#    for(i in 1:nrow(M)){A[M[i,1],M[i,2]]=A[M[i,2],M[i,1]]=1}
#    diag(A)=0
#    while(TRUE){
#      B=A
#      for(i in 1:n){
#        for(j in 1:n){
#          for(k in 1:n){
#            if(A[i,j]>A[i,k]+A[k,j]){
#              A[i,j]=A[i,k]+A[k,j]
#            }
#          }
#        }
#      }
#      if(identical(B,A)){break}else{B=A} ##check whether it is reasonable
#    }
#    return(A)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix Rcpp_shortest(NumericMatrix M){
#    M=M-1;
#    int n=max(M(_,1))+1;
#    NumericMatrix A(n,n);
#    A.fill(R_PosInf);
#    A.fill_diag(0);
#    for(int i=0;i<M.nrow();i++){
#      A(M(i,0),M(i,1))=1;
#      A(M(i,1),M(i,0))=1;
#    }
#    while(true){
#      NumericMatrix B=clone(A);
#      for(int i=0;i<n;i++){
#        for(int j=0;j<n;j++){
#          for(int k=0;k<n;k++){
#            if(A(i,j)>A(i,k)+A(k,j)){
#              A(i,j)=A(i,k)+A(k,j);
#            }
#          }
#        }
#      }
#      if(sum(A!=B)==0){break;}else{B=clone(A);}
#    }
#    return A;
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp22031)
data(data)
identical(R_shortest(data),Rcpp_shortest(data)) ##consider whether the returned results are consistent
library(rbenchmark)
benchmark(R_shortest(data),Rcpp_shortest(data),
          columns=c("test","replications","elapsed","relative")) #compare the speed

