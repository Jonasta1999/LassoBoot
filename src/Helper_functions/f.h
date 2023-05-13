double f(arma::vec v) {
  v = v.elem(find_finite(v));
  return sqrt(sum(square(v)) / std::max(1.0, v.n_elem - 1.0));
}