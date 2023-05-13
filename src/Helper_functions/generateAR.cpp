#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
#include<scale.h>
#include <generateX.h>


// [[Rcpp::export]]
arma::mat generateAR(int n, int p, double kappa, double rho){
  
  arma::mat y(n+1, p, fill::zeros);
  arma::mat e = generateX(n+1, p, kappa);
  
  for(int i = 1; i < (n+1); i++){
    y.row(i) = rho*y.row(i-1) + e.row(i);
  }
  return(y.rows(1,n));
}