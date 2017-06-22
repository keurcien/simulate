#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' Peak detection
//' 
//' \code{detect_cpp} returns a numeric vector determining where switches
//' are located.
//' 
//' @param stat a numeric vector.
//' @param tol an integer.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericVector detect_cpp(const NumericVector &stat, double threshold){
  int n = stat.size();  
  NumericVector bounds(n);
  bool switch_indicator = 0;
  for (int i = 0; i < n - 1; i++){
    if ((stat[i] > threshold) && (stat[i + 1] > threshold) && (switch_indicator == 0)){
      bounds[i] = 1;
      switch_indicator = 1;
    }
    if ((stat[i] > threshold) && (stat[i + 1] < threshold) && (switch_indicator == 1)){
      bounds[i] = 1;
      switch_indicator = 0;
    }
  }
  return(bounds);
}