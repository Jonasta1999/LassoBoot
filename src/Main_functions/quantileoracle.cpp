#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include<../Helper_functions/quantiles_arma.h>
#include<../Helper_functions/matrep.h>
#include<../Helper_functions/scale.h>
#include<../Helper_functions/generateX.h>
#include<../Helper_functions/generateAR.h>


// [[Rcpp::export]]
arma::vec quantileoracle(int n, int p, double kappa, double rho, const arma::vec &level, float sdev = 1,  int nbr_bootstraps=1000){

  arma::vec eff_noise(nbr_bootstraps, fill::zeros);

  for(int i = 0; i < nbr_bootstraps; i++){

    arma::vec noise = sdev * randn(n);
    arma::mat X = generateAR(n, p, kappa, rho);

    eff_noise(i) = 2 * max(abs(trans(X) * noise)) / X.n_rows;

  }
  arma::vec oracle_lambda(level.n_elem, fill::zeros);
  for(int j = 0; j < level.n_elem; j++){
    oracle_lambda(j) = arma::as_scalar(quantile_arma(eff_noise, level(j)));
  }

  return(oracle_lambda);
}





