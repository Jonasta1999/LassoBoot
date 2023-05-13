#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include<../Helper_functions/matrep.h>
#include<../Helper_functions/scale.h>


// [[Rcpp::export]]
arma::mat generateX(int n, int p, double kappa, bool scaling = false){
  
  arma::mat sigma = (1 - kappa) * arma::eye(p,p) +  arma::ones(p,p)*kappa;
  arma::colvec mu(p, fill::zeros);
  arma::mat X(n, 1, fill::zeros);
  arma::mat res = arma::mvnrnd(mu, sigma, n).t();
  if(scaling == true){
    res = scaleData(res, false) / ( sqrt(res.n_rows-1) * sqrt(res.n_rows) );
  }
  return(res);
}