#include <RcppEigen.h>
#include <Rcpp.h>


// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;
using namespace Rcpp;
//' Eigen_C
//' @param As A sysmetric matrix
//' @export
//[[Rcpp::export]]

List Eigen_C(NumericMatrix As){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(As));
  SelfAdjointEigenSolver<MatrixXd> es(A);
  return List::create(Named("values") = es.eigenvalues().reverse(),
                      Named("vectors") = es.eigenvectors().rowwise().reverse());

}

//' Eigen_C_value
//' @param As A sysmetric matrix
//' @export
//[[Rcpp::export]]
Eigen::VectorXd Eigen_C_value(NumericMatrix As){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(As));
  SelfAdjointEigenSolver<MatrixXd> es(A);
  return  es.eigenvalues().reverse();
}


/* calculate t(A)KA*/
  /* this function computes everything in through elementwise style*/
  /* seems very slow in Rcpp syxtax */
  /* works well in C*/

