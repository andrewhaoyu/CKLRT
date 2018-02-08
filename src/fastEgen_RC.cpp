#include <RcppEigen.h>
#include <Rcpp.h>


// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;
typedef  Map<VectorXd>  MapVecd;
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

//' MatMult_C
//' @param A first matrix
//' @param B second matrix
//' @export
// [[Rcpp::export]]
SEXP MatMult_C(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}


//' Sum_C
//' @param AA Vector
//' @export
// [[Rcpp::export]]
double Sum_C(NumericVector AA){
  const MapVecd A(as<MapVecd>(AA));
  double result= A.sum();
  return result;
}

//' ColSum_C
//' @param AA Matrix
//' @export
// [[Rcpp::export]]
NumericVector ColSum_C(NumericMatrix AA){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  NumericVector result;
  result = A.colwise().sum();
  return result;
}

//' MatrixRowMax_C
//' @param AA Matrix
//' @export
// [[Rcpp::export]]
NumericVector MatrixRowMax_C(NumericMatrix AA){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  NumericVector result;
  result = A.rowwise().maxCoeff();
  return result;
}



//' Elementwisesquare_C
//' @param AA Matrix
//' @export
// [[Rcpp::export]]
NumericVector Elementwisesquare_C(NumericMatrix AA){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  NumericVector result;
  result = A.array().square();
  return result;
}




//' VecMultMat_C
//' @param A Vector
//' @param B Matrix
//' @export
// [[Rcpp::export]]
SEXP VecMultMat_C(Eigen::VectorXd A,Eigen::MatrixXd  B){
  Eigen::VectorXd C = A.transpose()*B;
  return Rcpp::wrap(C);
}




//' ColSumtwomatrix_C
//' @param AA Matrix
//' @param BB Matrix
//' @export
// [[Rcpp::export]]
NumericVector ColSumtwomatrix_C(NumericMatrix AA,NumericMatrix BB){
  NumericVector result;
  result = ColSum_C(AA)+ColSum_C(BB);
  return result;
}





//' MatrixPlus_C
//' @param A First Matrix
//' @param B Second Matrix
//' @export
// [[Rcpp::export]]
SEXP MatrixPlus_C(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A + B;
  return Rcpp::wrap(C);
}


//' NumxMatrix_C
//' @param A Number
//' @param B Matrix
//' @export
// [[Rcpp::export]]
SEXP NumxMatrix_C(double A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

/* calculate t(A)KA*/
  /* this function computes everything in through elementwise style*/
  /* seems very slow in Rcpp syxtax */
  /* works well in C*/

