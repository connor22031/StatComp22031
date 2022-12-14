% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp22031R.R
\name{benchmarks}
\alias{benchmarks}
\title{Benchmark R and Rcpp functions.}
\description{
Use R package \code{rbenchmark} and \code{microbenchmark} to compare the performance of R functions (\code{R_shortest}) and Cpp functions (\code{Rcpp_shortest}).
}
\examples{
\dontrun{
data(data)
attach(data)
benchmark(R_shortest(data),Rcpp_shortest(data),
columns=c("test","replications","elapsed","relative"))
microbenchmark::microbenchmark(gibbs_cpp(5000,0,0,1,1,0.9),
gibbs_r(5000,0,0,1,1,0.9))
}
}
