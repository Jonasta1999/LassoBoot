
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// main function
//'  Calculate the estimated quantiles of the effective noise
//'
//' 'quantileestimated' lets you calculate the estimated quantiles using the criterion quantiles.
//'
//' @param Y Response variable vector.
//' @param X Covariates matrix.
//' @param level Vector of levels.
//' @param estimator_grid Lasso estimator grid matrix.
//' @param nbr_bootstraps Number of bootstrap simulations. Defaults to 1000.
//' @param criterion_quantiles The criterion quantiles computed using the 'criterionquantiles' function.
//' @return quantilesestimated Returns the value of the estimated quantile of the Lasso's effective noise.
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
//'                                   bs_type = "stationary_block")
//'quantile_e <- quantileestimated(Y = Y,
//'                                 X = X,
//'                                 level = 0.95,
//'                                 estimator_grid,
//'                                 criterionquantiles = criterion_q)
//'
//'
// [[Rcpp::export]]
arma::vec quantileestimated(const arma::mat &Y, const arma::vec &X, const arma::vec &level, arma::vec &lambda_grid, arma::mat &criterion_quantiles){
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
