#include <Rcpp.h>
using namespace Rcpp;

//' @title Shortest path model between cities
//' @description Output the connectivity between each city
//' @param M a direct connection situation matrix of size \code{n*2}
//' @return a matrix of size \code{n*n}
//' @examples
//' \dontrun{
//' data(data)
//' attach(data)
//' R_shortest(data)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix Rcpp_shortest(NumericMatrix M){
  M=M-1;
  int n=max(M(_,1))+1; 
  NumericMatrix A(n,n);
  A.fill(R_PosInf);
  A.fill_diag(0);
  for(int i=0;i<M.nrow();i++){
    A(M(i,0),M(i,1))=1;
    A(M(i,1),M(i,0))=1;
  }
  while(true){
    NumericMatrix B=clone(A);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        for(int k=0;k<n;k++){
          if(A(i,j)>A(i,k)+A(k,j)){
            A(i,j)=A(i,k)+A(k,j);
          }
        }
      }
    }
    if(sum(A!=B)==0){break;}else{B=clone(A);}
  }
  return A;
}

//' @title A binary normal Gibbs sampler using Rcpp
//' @description Generate a binary normal sample with Gibbs
//' @param N the sample size
//' @param mu1 the mean value of one of the sample
//' @param mu2 the mean value of another sample
//' @param sigma1 the variance of one of the sample
//' @param sigma2 the variance of another sample
//' @param rho the sample correlation coefficient
//' @return A binary normal random sample of size \code{N}
//' @examples
//' \dontrun{
//' gibbs_cpp(5000,0,0,1,1,0.9)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbs_cpp(int N,double mu1,double mu2,double sigma1,double sigma2,double rho){
  NumericMatrix X(N,2);
  double s1,s2;
  s1=sqrt(1-pow(rho,2))*sigma1;
  s2=sqrt(1-pow(rho,2))*sigma2;
  X(0,0)=mu1;
  X(0,1)=mu2;
  double x1,x2,m1,m2;
  for(int i=1;i<N;i++){
    x2=X(i-1,1);
    m1=mu1+rho*(x2-mu2)*sigma1/sigma2;
    X(i,0)=rnorm(1,m1,s1)[0];
    x1=X(i,0);
    m2=mu2+rho*(x1-mu1)*sigma2/sigma1;
    X(i,1)=rnorm(1,m2,s2)[0];
  }
  return(X);
}