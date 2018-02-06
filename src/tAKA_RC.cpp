#include <Rcpp.h>


using Rcpp::NumericMatrix;
//[[Rcpp::export]]

Rcpp::NumericMatrix tAKA_RC(NumericMatrix A,NumericMatrix K){
  double sum = 0.0;

  int Anc = A.ncol();
  int Anr = A.nrow();
  NumericMatrix ret(Anc,Anc);

  for(int i=0;i<Anc;i++){
    for(int j=0;j<(i+1);j++){
      /*the ith row and jth column of the ret*/
        /*the ith column X transpose times K times the jth column of the X */
        sum = 0.0;
        /* One vector times K times one Vector*/
          for(int k=0;k<Anr;k++){
            for(int l=0;l<Anr;l++){
              sum += A(k,i)*K(k,l)*A(l,j);
            }
          }
        ret(i,j) = sum;
        /* ret is sysmetric */
          ret(j,i) = sum;

    }
  }
  return ret;
}

