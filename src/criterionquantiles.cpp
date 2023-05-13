#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <regex>

using namespace Rcpp;
using namespace arma;

#include"Helper_functions/quantiles_arma.h"
#include"Helper_functions/matrep.h"
#include"block_bootstrap/stationary_bootstrap.h"
#include"block_bootstrap/block_bootstrap.h"

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
//'  Calculate criterion quantiles
//'
//' 'criterionquantiles' lets you calculate the criterion quantiles using a bootstrap.
//'
//' @param y Response variable vector.
//' @param X Covariates matrix.
//' @param level Vector of levels.
//' @param estimator_grid Lasso estimator grid matrix.
//' @param nbr_bootstraps Number of bootstrap simulations. Defaults to 1000.
//' @param w Number of blocks. If bs_type = "multiplier" then w is not applicable. If bs_type = "stationary", then w is the mean of number of blocks used for each bootstrap sample. Defaults to 10.
//' @param bs_type Type of bootstrap technique to be used. The following methods are available for use: "multiplier", "stationary_block" and "moving_block". Any specified character length is sufficient. Defaults to "multiplier".
//' @return criterion_quantiles Returns a matrix of quantile at each level with each row being the different values from the grid of lambdas.
//'
//'@examples
//'lasso <- lasso(Y = Y,
//'                X = X,
//'                out = 1,
//'                numlam = 100,
//'                maxIter = 1e10,
//'                opTol = 1e-10,
//'                standardize = 0)
//'
//'#Remove equidistant lambda row
//'BETA <- lasso %>% head(-1)
//'
//'criterion_q <- criterionquantiles(Y = Y,
//'                                   X = X,
//'                                   level = c(0.90, 0.95, 0.99),
//'                                   estimator_grid = BETA,
//'                                   nbr_bootstraps = 10000,
//'                                   w = 20,
//'                                   bs_type = "det")
// [[Rcpp::export]]
arma::mat criterionquantiles(const arma::vec& y, const arma::mat& X,
                             const arma::vec& level, const arma::mat& estimator_grid,
                             int nbr_bootstraps = 1000, int w = 10, std::string bs_type = "multiplier") {

  arma::mat criterion_quantiles(estimator_grid.n_cols, level.n_elem, fill::zeros);
  arma::mat residual_ext = matrep(y, estimator_grid.n_cols, X.n_rows, estimator_grid.n_cols) - (X * estimator_grid);
  arma::mat residual_large(residual_ext.n_rows, residual_ext.n_cols* nbr_bootstraps, fill::zeros);

  bs_type.insert(0, "^");
  if(std::regex_search("stationary_block", std::regex(bs_type)) ){
    residual_large = stationary_bootstrap(residual_ext, nbr_bootstraps, w = w);
  }
  else if(std::regex_search("moving_block", std::regex(bs_type))){
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
    throw std::invalid_argument("Block bootstrap type must be either 'multiplier', 'stationary_block' or 'moving_block.");
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

