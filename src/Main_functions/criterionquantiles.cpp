#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <regex>

using namespace Rcpp;
using namespace arma;

#include<../Helper_functions/quantiles_arma.h>
#include<../Helper_functions/matrep.h>
#include<../../block_bootstrap/stationary_bootstrap.h>
#include<../../block_bootstrap/block_bootstrap.h>

// [[Rcpp::depends(RcppArmadillo)]]


// function to find column maximums
arma::vec colMaxs(const arma::mat& X) {
  int n = X.n_cols;
  arma::vec out(n);
  for (int i = 0; i < n; i++) {
    out(i) = X.col(i).max();
  }
  return out;
}

// main function
// [[Rcpp::export]]
arma::mat criterionquantiles(const arma::vec& y, const arma::mat& X, 
                             const arma::vec& level, const arma::mat& estimator_grid, 
                             int nbr_bootstraps = 1000, int w = 10, std::string bs_type = "probabilistic") {
  int n = y.n_elem;
  int p = X.n_cols;
  int l = level.n_elem;
  int e = estimator_grid.n_cols;
  int X_rows = X.n_rows;
  
  arma::mat criterion_quantiles(estimator_grid.n_cols, level.n_elem, fill::zeros);
  arma::mat residual_ext = matrep(y, estimator_grid.n_cols, X.n_rows, estimator_grid.n_cols) - (X * estimator_grid);
  arma::mat residual_large(residual_ext.n_rows, residual_ext.n_cols* nbr_bootstraps, fill::zeros);
  
  bs_type.insert(0, "^");
  if(std::regex_search("probabilistic", std::regex(bs_type)) ){
    residual_large = stationary_bootstrap(residual_ext, nbr_bootstraps, w = w);
  }
  else if(std::regex_search("deterministic", std::regex(bs_type))){
    residual_large = block_bootstrap(residual_ext, nbr_bootstraps, w = w);
  }
  else if(std::regex_search("multiplier", std::regex(bs_type))){
    residual_large = repmat(residual_ext, 1, nbr_bootstraps);
    residual_large = residual_large % arma::kron(reshape(randn(X.n_rows * nbr_bootstraps), X.n_rows, nbr_bootstraps), reshape(ones(1, estimator_grid.n_cols), 1, estimator_grid.n_cols));
    
    arma::mat effective_noise_matrix = trans(X) * residual_large;
    arma::vec effective_noise = 2 * colMaxs(abs(effective_noise_matrix)) / X.n_rows;
    
    for (int i = 0; i < estimator_grid.n_cols; i++) {
      arma::uvec indexes = arma::regspace<arma::uvec>(i*nbr_bootstraps, i*nbr_bootstraps + nbr_bootstraps-1);
      for (int j = 0; j < level.n_elem; j++) {
        criterion_quantiles(i, j) = quantile_arma(effective_noise.rows(indexes), level(j));
      }
    }
    return(criterion_quantiles);
    
  }
  else{
    throw std::invalid_argument("Block bootstrap type must be either 'probabilistic', 'deterministic' or 'multiplier.");
  }
  
  
  
  arma::mat effective_noise_matrix = trans(X) * residual_large;
  arma::vec effective_noise = 2 * colMaxs(abs(effective_noise_matrix)) / X.n_rows;
  
  for (int i = 0; i < estimator_grid.n_cols; i++) {
    arma::uvec indexes = arma::regspace<arma::uvec>(i*nbr_bootstraps, i*nbr_bootstraps + nbr_bootstraps-1);
    for (int j = 0; j < level.n_elem; j++) {
      criterion_quantiles(i, j) = quantile_arma(effective_noise.rows(indexes), level(j));
    }
  }
  return(criterion_quantiles);
}

