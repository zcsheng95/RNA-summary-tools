#include <RcppEigen.h>

inline Eigen::MatrixXd AtA(const Eigen::MatrixXd &A) {
  int n(A.cols());
  return Eigen::MatrixXd(n, n)
      .setZero()
      .selfadjointView<Eigen::Lower>()
      .rankUpdate(A.adjoint());
}
