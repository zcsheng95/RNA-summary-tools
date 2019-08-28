#include "util.cpp"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

//' myLSrcpp
//' The function illustrates how to use RcppEigen to fit the model
//' Y = X beta + epsilon
//' What is returned is the vector betahat=LSE(beta)
//' @author Kouros Owzar
//' @param Y the vector of response
//' @param X the matrix of covariates (does not include intercept)
//' @return \code{betahat}
//' @export
// [[Rcpp::export]]
Eigen::VectorXd myLSrcpp(const Eigen::VectorXd Y, const Eigen::MatrixXd X) {

  /* This function returns the LSE vector */

  // Get number of observations
  int n = X.rows();
  // Get number of variables
  int m = X.cols();

  // Create Intercept vector
  Eigen::VectorXd onevec(n);
  onevec.fill(1);

  // Create design matrix by augmenting the intercept vector
  // to the covariate matrix
  Eigen::MatrixXd X0(n, m + 1);
  X0 << onevec, X;

  // Excercise (what happens here)?
  const Eigen::LLT<Eigen::MatrixXd> llt(AtA(X0));
  // Calculate the LSE vector
  const Eigen::VectorXd betahat(llt.solve(X0.adjoint() * Y));

  return (betahat);
}
