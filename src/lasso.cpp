#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// main function
//'  Compute lasso
//'
//' 'lasso' lets you compute the lasso estimator using the coordinate descent algorithm.
//'
//' @param Y Response variable vector.
//' @param X Covariates matrix.
//' @param out 1: Outputs betas. 2: Outputs residuals.
//' @param numlam Number of lambdas in equidistant lambda grid.
//' @param maxIter Maximum number of iterations.
//' @param optTol The optimization tolerance as an indicatin for convergence.
//' @param standardize Whether or not data should be standardized.
//' @return lasso Returns lasso estimates for each lambda.
//'
//'@examples
//'BETA <- lasso(X = X, Y = Y, out = 1)
//'
// [[Rcpp::export]]
arma::mat lasso(const arma::mat& X, const arma::vec& Y, int out, int numlam, double maxIter, double optTol, double standardize) {

  int n = X.n_rows, p = X.n_cols;
  double m;
  arma::mat X1 = X;
  arma::vec Y1 = Y;

  if (standardize == 1){
    arma::mat M1 = mean(X,0);
    arma::mat S1 = stddev(X,0,0);
    X1 = (X-repmat(M1, n, 1))/repmat(S1, n, 1);
    double M2 = mean(Y);
    Y1 = Y - M2;
  } else {
    X1 = X;
    Y1 = Y;
  }
  arma::mat XX = X1.t()*X1;
  arma::mat Xy = X1.t()*Y1;
  arma::vec beta(p, fill::zeros);

  arma::mat XX2 = XX*2/n;
  arma::mat Xy2 = Xy*2/n;

  arma::mat BETA(p,numlam);
  arma::mat RES(n,numlam);
  // definition of lambda equidistant grid
  double lam_max = 2*max(abs(X1.t()*Y1)) / n;
  arma::vec Lambda1 = linspace(lam_max, 0.0, numlam+1);
  for (int l1 = 0; l1 < numlam; l1++){ // loop through Lambda1
    m = 0;
    double lambda1 = Lambda1[l1];
    while (m < maxIter){  // while loop until convergence beta^(i)
      arma::vec betaold = beta;
      for(int j = 0; j < p; j++){ // loop through beta^(i)_j
        // compute the Shoot and Update the variable
        double S0 = sum(XX2.row(j) * beta) - XX2(j,j) * beta[j] - Xy2[j];
        if (S0 > lambda1){
          beta[j] = (lambda1 - S0)/XX2(j,j);
        }
        if (S0 < -lambda1){
          beta[j] = (-lambda1 - S0)/XX2(j,j);
        }
        if (fabs(S0) <= lambda1) {
          beta[j] = 0;
        }
      } // end loop beta^(i)_j
      m = m + 1;
      if (sum(abs(beta - betaold)) < optTol) {
        m = maxIter + 1;
      }
    } // end while loop
    BETA.col(l1) = beta;//BETA  store beta
    if (out == 2){
      RES.col(l1) = Y1 - X1*beta;
    }
  }  // end loop Lambda1
  arma::mat OUT;
  if (out == 1){
    arma::vec rev_lam = reverse(Lambda1);
    rev_lam.shed_row(0);
    OUT = join_cols(reverse(BETA, 1), rev_lam.t());
  }
  if (out == 2){
    //OUT = RES;
    arma::vec rev_lam = reverse(Lambda1);
    rev_lam.shed_row(0);
    OUT = join_cols(RES, rev_lam.t());
    //OUT = rev_lam;
  }
  return OUT;
}
