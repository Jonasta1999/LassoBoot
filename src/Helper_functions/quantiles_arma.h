#ifndef QUANTILES_ARMA_H
#define QUANTILES_ARMA_H

inline double quantile_arma(arma::vec vec, double p)
{
  int n = vec.size();
  arma::uvec sorted_index = arma::sort_index(vec);

  // Index of quantile to be returned
  int index = (int)ceil(p * n) - 1;

  // Return the quantile
  return vec[sorted_index[index]];
}

#endif // QUANTILES_ARMA_H
