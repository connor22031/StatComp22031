% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp22031R.R
\name{bootstrap}
\alias{bootstrap}
\title{The bootstrap method}
\description{
Use R package \code{boot} and \code{bootstrap} to implement bootstrap method
}
\examples{
\dontrun{
t=c(3,5,7,18,43,85,91,98,100,130,230,487)
p=function(x,i){
 1/mean(x[i])
}
boot(t,statistic = p,R=2000)
}
}
