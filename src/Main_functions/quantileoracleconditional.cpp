
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include<../Helper_functions/matrep.h>
#include<../Helper_functions/colmax.h>
#include<../Helper_functions/quantiles_arma.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec quantileoracleconditional(const arma::mat &X, const arma::vec &level, float sdev = 1,int nbr_bootstraps = 1000){
  
  arma::mat noise = randn(X.n_rows, nbr_bootstraps) * sdev;
  arma::mat effective_noise_matrix = trans(X) * noise;
  arma::vec effective_noise = 2 * colMaxs(abs(effective_noise_matrix)) / X.n_rows;
  arma::vec oracle_lambda(level.n_elem, fill::zeros);
  
  for(int i = 0; i < level.n_elem; i++){
    oracle_lambda(i) = arma::as_scalar(quantile_arma(effective_noise, level(i)));
  }
  
  return(oracle_lambda);

}
  