
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

arma::vec quantileestimated(const arma::mat &Y, const arma::vec &X, const arma::vec &level, arma::vec &lambda_grid, arma::mat &criterion_quantiles, int nbr_lambda = 100, int nbr_bootstraps = 1000){
  arma::vec selected_lambda(level.n_elem, arma::fill::zeros);
  
  for(int ll = 0; ll < level.n_elem; ll++){
    
    int selected_index = lambda_grid.n_elem - 1;
    while ((selected_index > 0) && (criterion_quantiles(selected_index, ll) <= lambda_grid(selected_index)))
      selected_index--;
    if(selected_index > 0)
      selected_index = std::min(selected_index + 1, (int)lambda_grid.n_elem - 1);
    if((selected_index == 0) && (criterion_quantiles(0, ll) > lambda_grid(0)))
      selected_index = 1;
    selected_lambda(ll) = criterion_quantiles(selected_index, ll);
  }
  return selected_lambda;
}