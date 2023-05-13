#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat matrep(const arma::mat& Y, const int& k, const int& n, const int& m){
  arma::mat Y_rep = arma::repmat(Y,1,k);
  return(arma::reshape(Y_rep,n,m));
}






