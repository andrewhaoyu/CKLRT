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
NumericVector VecMultMat_C(Eigen::VectorXd A,Eigen::MatrixXd  B){
  Eigen::VectorXd C = A.transpose()*B;
  return Rcpp::wrap(C);
}

//' Vecplus_C
//' @param A Vector
//' @param B Vector
//' @export
// [[Rcpp::export]]
NumericVector Vecplus_C(Eigen::VectorXd A,Eigen::VectorXd B){
  Eigen::VectorXd C = A+B;
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

//' ifelsetest_C
//' @param x Vector
//' @export
//[[Rcpp::export]]
NumericVector ifelsetest_C( NumericVector x){
  return Rcpp::wrap( ifelse( x < 0, 0, x ));
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


//' LR0_fixRho_C
//' @param LamdasR Lamda Number
//' @param muR mu vector
//' @param w1R w1 vector
//' @param w2R w2 vector
//' @param nminuspx n-px
//' @export
// [[Rcpp::export]]
/*NumericMatrix LR0_fixRho_C(NumericVector LamdasR,
             NumericVector muR,
             NumericMatrix w1R,
             NumericVector w2R,
             int nminuspx)*/
NumericVector LR0_fixRho_C(NumericVector LamdasR,
                           NumericVector muR,
                           NumericMatrix w1R,
                           NumericVector w2R,
                           int nminuspx)
  {
  const Map<MatrixXd> w1(as<Map<MatrixXd> >(w1R));
  const Map<VectorXd> w2(as<Map<VectorXd> >(w2R));
  const Map<VectorXd> mu(as<Map<VectorXd> >(muR));
  const Map<VectorXd> Lamdas(as<Map<VectorXd> >(LamdasR));
  int length_lamda = Lamdas.size();
  int N = w2.size();
  double lam;
  Eigen::VectorXd lammu_con;
  Eigen::VectorXd lammu_case;
  Eigen::VectorXd Dn;
  Eigen::VectorXd Nn;
  NumericVector temp;
  NumericMatrix result(N,length_lamda);
  int i = 0;
  /*for(int i=0;i<length_lamda;i++){*/
    lam = Lamdas[i];
    lammu_con = 1/(lam*mu).array();
    lammu_case = 1- lammu_con.array();
    Dn = (lammu_con).transpose()*w1+w2;
    Nn = lammu_case*w1;
    temp = nminuspx*(1+Nn.array()/Dn.array()).log()-
      (1+lam*mu.array()).log().sum();
    /*result(_,i) = ifelsetest_C(temp);*/
  /*}*/
  return Rcpp::wrap(ifelsetest_C(temp));
}

/* calculate t(A)KA*/
  /* this function computes everything in through elementwise style*/
  /* seems very slow in Rcpp syxtax */
  /* works well in C*/

