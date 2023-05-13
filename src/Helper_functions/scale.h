
arma::vec colMeans(const arma::mat & x) {
  arma::vec means = mean(x, 0);
  return means;
}


double f(const arma::vec & v) {
  arma::vec v_non_na = v.elem(find_finite(v));
  arma::vec m(2, fill::zeros);
  m(0) = 1;
  m(1) = v_non_na.n_elem - 1.0;
  return sqrt(sum(square(v_non_na)) / max(m));
}


arma::mat scaleData(const arma::mat & x, const bool center = true, const bool scale = true) {
  arma::mat x_scaled = x;
  int nc = x.n_cols;
  if (center) {
    arma::vec center = colMeans(x_scaled);
    
    x_scaled.each_col() -= center;
  }
  if (scale) {
    arma::vec scale(nc);
    for (int j = 0; j < nc; j++) {
      scale(j) = f(x_scaled.col(j));
    }
    
    std::cout << x_scaled <<"\n" << trans(scale)<<"\n";
    
    for(int i = 0; i < x_scaled.n_rows; i++){
      for(int k = 0; k < x_scaled.n_cols; k++){
        x_scaled(i,k) /= scale(k); 
      }
    }
  }
  return x_scaled;
}









