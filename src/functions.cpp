#include <Rcpp.h>
using namespace Rcpp;

// Something broken in the devtools config without some c++ function here
// // [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}