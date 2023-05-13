
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include"Helper_functions/matrep.h"
#include"Helper_functions/colmax.h"
#include"Helper_functions/quantiles_arma.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// main function
//'  Calculates the oracle quantile
//'
//' 'quantileoracleconditinal' lets you calculate the oracle quantiles dependent on X.
//'
//' @param X Covariates matrix.
//' @param level Vector of levels.
//' @param sdev The specified standard deviation for the noise.
//' @param nbr_bootstraps Number of bootstrap simulations. Defaults to 1000.
//' @return quantileoracleconditional Returns the value of the approximated quantile of the lasso's effective noise dependent on X.
//'
//'@examples
//'
//'quantileoracleconditional(X = X,
//'                           level = c(0.95, 0.99),
//'                           sdev = 1)
//'
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
