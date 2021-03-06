#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' Jumps
//' 
//' \code{jumps_from_map} returns a numeric vector determining where jumps
//' are located.
//' 
//' @param map a numeric vector.
//' @param lambda a numerical value specifying the parameter of the exponential
//' law.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericVector jumps_from_map(const NumericVector &map, const double &lambda){
  int n = map.size();  
  NumericVector jumps(n);
  double g = 0;
  double p = 0;
  double ber = 0;
  for (int i = 1; i < n; i++){
    g += map[i] - map[i - 1];
    p = R::pexp(g, lambda, 1, 0); // g in Morgans
    ber = R::rbinom(1.0, p);
    if (ber == 1){
      g = 0;
      jumps[i] = 1;
    }
  }
  return(jumps);
}

//' Chunks
//' 
//' \code{ancestry_chunks} converts jumps to chunks.
//' 
//' @param jumps a numeric vector obtained with \code{jumps_from_map}.
//' 
//' @return The returned value is a logical vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
LogicalVector ancestry_chunks(const NumericVector &jumps){
  int n = jumps.size();  
  LogicalVector chunks(n);
  bool fbool = true;
  for (int i = 0; i < n; i++){
    if (jumps[i] == 1){
      fbool = !fbool;
    }
    chunks[i] = fbool;
  }
  return(chunks);
}

//' Convert haplotypes to genotypes
//' 
//' \code{haplo_to_geno} converts haplotypes to genotypes.
//' 
//' @param H a numerical matrix.
//' 
//' @return The returned value is a numerical matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericMatrix haplo_to_geno(const NumericMatrix &H){
  int nSNP = H.nrow();
  int nIND = H.ncol() / 2;
  NumericMatrix G(nSNP, nIND);
  for (int j = 0; j < nIND; j++){
    for (int i = 0; i < nSNP; i++){
      G(i, j) = H(i, 2 * j) + H(i, 2 * j + 1); 
    }
  }
  return G;
}
