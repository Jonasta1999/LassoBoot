arma::mat stationary_bootstrap(const arma::mat data, int B = 1000, int w = 10) {
  int t = data.n_rows;
  int k = data.n_cols;
  
  // Define the probability of a new block
  double p = 1.0/w;
  
  // Set up bootstrap data and indices
  arma::mat indices(t, B, fill::zeros);
  
  
  // Initial positions
  indices.row(0) = arma::ceil(t * randu<rowvec>(B));
  
  // Set up the random numbers
  arma::mat select(t, B, fill::randu);
  select.transform([p](double val) { return val < p ? 1.0 : 0.0; });
  arma::umat select_indices = arma::find(select);
  indices(select_indices) = ceil( randu( select_indices.n_elem ) * t );
  
  for(int i = 1; i < indices.n_rows; i++){
    for(int j = 0; j < indices.n_cols; j++){
      if(indices(i,j) == 0){
        indices(i,j) = indices(i-1,j) + 1.0;
      }
    }
  }
  indices.transform([t](double val) { return val > t ? val - t: val; });
  indices.reshape(t*B, 1);
  
  // Create new bsdata matrix identical to indices
  arma::mat bsdata(t*B, k, fill::zeros);
  
  // Convert all elements to the corresponding element from 'data' based on indices
  for(int i = 0; i < t*B; i++){
    bsdata.row(i) = data.row(indices(i,0)-1);
  }
  bsdata.reshape(t, k*B);
  return(bsdata);
}