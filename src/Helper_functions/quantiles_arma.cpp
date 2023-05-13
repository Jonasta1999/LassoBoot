#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
double quantile_arma(const vec& x, double p) {
  vec sorted_x = sort(x);
  int n = x.n_elem;
  double index = p * (n - 1);
  int lower_idx = std::floor(index);
  int upper_idx = std::ceil(index);
  
  if (lower_idx == upper_idx) {
    return sorted_x(lower_idx);
  } else {
    double lower_weight = 1 - (index - lower_idx);
    double upper_weight = 1 - (upper_idx - index);
    return lower_weight * sorted_x(lower_idx) + upper_weight * sorted_x(upper_idx);
  }
}
