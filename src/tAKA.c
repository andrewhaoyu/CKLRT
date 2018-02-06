#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


#define ALMOST_ZERO 1e-16
#define NUM_ZERO 1e-100
#define ERROR_SINGULAR_MATRIX 1
#define CHECK_MEM(obj) if (obj == NULL) {Rprintf("ERROR: allocating memory \n"); error("1");}

static void print_dVec(vec, n, name)
  double *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %g ", vec[i]);
  }
  printf("\n \n");
}
static void print_iVec(vec, n, name)
  int *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %d ", vec[i]);
  }
  printf("\n \n");
}


static void print_dMat(mat, nr, nc, name)
  double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  printf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) printf(" %g ", mat[i][j]);
    printf("\n");
  }
  printf("\n \n");
}

/* Function to allocate memory for a double vector */
static double * dVec_alloc(n, initFlag, initVal)
  int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */

/* Function to allocate a double matrix */
static double ** dMat_alloc(nrow, ncol, initFlag, initVal)
  int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: dMat_alloc */

/* Function to free a matrix */
static void matrix_free(x, n)
  void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

/* Function to check for non-finite values */
static int all_finite(vec, n)
  double *vec;
int n;
{
  int i;

  for (i=0; i<n; i++) {
    if (!R_FINITE(vec[i])) return(0);
  }

  return(1);

} /* END: all_finite */

/* Multiply to matrices (matrix mult)  */
static void matrixMult(m1, m1_nr, m1_nc, m2, m2_nc, ret)
  double **m1, **m2, **ret;
int m1_nr, m1_nc, m2_nc;
{
  int i, j, k;
  double sum;

  for (i=0; i<m1_nr; i++) {
    for (j=0; j<m2_nc; j++) {
      sum = 0.;
      for (k=0; k<m1_nc; k++) sum += m1[i][k]*m2[k][j];
      ret[i][j] = sum;
    }
  }

} /* END: matrixMult */


/* two matrix minues  */
static void matrixminus(double **mat1,double**mat2,int nr,int nc,double **ret){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      ret[i][j] = mat1[i][j]-mat2[i][j];
    }
  }
}/* END: matrixminus */





/* Function to compute Xy for matrix X and vector */
static void X_y(X, nr, nc, y, ret)
  double **X, *y, *ret;
int nr, nc;
{
  int i, j;
  double sum, *p, *pret, *px;

  for (i=0, pret=ret; i<nr; i++, pret++) {
    sum  = 0.0;
    for (j=0, p=y, px=X[i]; j<nc; j++, p++, px++) {
      sum += *px * *p;
    }
    *pret = sum;
  }

} /* END: X_y */

/* Function for dot product of two vectors */
static double dotProd(v1, v2, n)
  double *v1, *v2;
int n;
{
  int i;
  double sum=0.0;

  for (i=0; i<n; i++) sum += v1[i]*v2[i];

  return(sum);

} /* END: dotProd */

static void fillMat(vec, nr, nc, addInt, out)
  double *vec, **out;
int nr, nc, addInt;
{
  int i, j, col=0, ii;

  if (addInt) {
    /* Intercept for first column */
    for (i=0; i<nr; i++) out[i][0] = 1.0;
    col = 1;
  }

  ii = 0;
  for (j=0; j<nc; j++) {
    for (i=0; i<nr; i++) {
      out[i][col] = vec[ii];
      ii++;
    }
    col++;
  }

} /* END: fillMat */


    /* fill the info matrix to the result*/
    static void fill_Info(Info,ret_info,Nparm)
  double **Info,*ret_info;
int Nparm;
{
  int i,j;
  for(j=0;j<Nparm;j++){
    for(i=0;i< Nparm;i++){
      ret_info[i*Nparm+j] = Info[i][j];
    }
  }


}

/*end fill_Info*/


/* Function for quadractic computation X^tkX */
static void QuadXKX(double **X,double ** K, int Xnr,int Xnc, double** ret){
  double sum = 0.0;

  for(int i=0;i<Xnc;i++){
    for(int j=0;j<(i+1);j++){
      /*the ith row and jth column of the ret*/
      /*the ith column X transpose times K times the jth column of the X */
      sum = 0.0;
      /* One vector times K times one Vector*/
      for(int k=0;k<Xnr;k++){
        for(int l=0;l<Xnr;l++){

          sum += X[k][i]*K[k][l]*X[l][j];
        }
      }
      ret[i][j] = sum;
      /* ret is sysmetric */
      ret[j][i] = sum;

    }
  }
}

void tAKA(double * AV, double* KV, int *pArow,int *pAcol, double * outV){
  double **A;
  double ** K;
  double ** out;
  int Arow = *pArow;
  int Acol = *pAcol;
  /*Rprintf("Allocate memory\n");*/
  A = dMat_alloc(Arow,Acol,0,0.0);
  K = dMat_alloc(Arow,Arow,0,0.0);
  out = dMat_alloc(Acol,Acol,0,0.0);
  fillMat(AV,Arow,Acol,0,A);
  /*print_dMat(A, Arow, Acol, "A");*/
  fillMat(KV,Arow,Arow,0,K);
  /*print_dMat(K, Arow, Arow, "K");*/
  /*Rprintf("Compute\n");*/
  QuadXKX(A,K,Arow,Acol,out);
  /*print_dMat(out, Acol, Acol, "out");*/
  fill_Info(out,outV,Acol);
  matrix_free((void **)A, Arow);
  matrix_free((void **)K, Arow);
  matrix_free((void **)out,Acol);
  return;
}
