arma::mat block_bootstrap(const arma::mat data, int B = 1000, int w = 10) {
  int t = data.n_rows;
  int k = data.n_cols;
  
  // Compute the number of blocks nedded
  double s = std::ceil(t/(double)w);
  
  // Generate the starting points
  arma::mat bs(s, B, fill::randu);
  bs.transform([t](double val) { return floor(val*t) + 1; });
  arma::mat indices(s*w, B, fill::zeros);
  int index = 0;
  
  // Adder is a variable that needs to be added each loop
  arma::mat adder = arma::repmat(arma::regspace(0, w-1), B, 1);
  adder = arma::reshape(adder, w, B);
  
  for(int i = 0; i < t; i += w){
    arma::vec rows_indices = arma::regspace(i, (i+w-1));
    arma::uvec rows_indices2 = arma::conv_to<arma::uvec>::from(rows_indices);
    
    arma::mat temp = arma::repmat(bs.row(index), w, 1);
    temp = temp + adder;
    
    
    indices.rows(rows_indices2) = temp;
    
    index++;
  }
  indices = indices.rows(0,t-1);
  indices.transform([t](double val) { return val > t ? val - t: val; });
  indices.reshape(t*B, 1);
  
  
  arma::mat bsdata(t*B, k, fill::zeros);
  for(int i = 0; i < t*B; i++){
    bsdata.row(i) = data.row(indices(i,0)-1);
  }
  bsdata.reshape(t, k*B);
  return(bsdata);
  
  
}