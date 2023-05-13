#ifndef COLMAX_H
#define COLMAX_H

inline arma::vec colMaxs(const arma::mat& X) {
  int n = X.n_cols;
  arma::vec out(n);
  for (int i = 0; i < n; i++) {
    out(i) = X.col(i).max();
  }
  return out;
}

#endif // COLMAX_H
