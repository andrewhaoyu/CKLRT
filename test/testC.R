A.row <- 3
A.col <- 3
A <- matrix(rnorm(A.row*A.col),A.row,A.col)
K <- matrix(rnorm(A.row*A.row),A.row,A.row)

#######test the MatMult_C function
result <- t(A)%*%K%*%A
result2 <- MatMult_C(MatMult_C(t(A),K),A)
library(microbenchmark)
microbenchmark(
  t(A)%*%A,
  crossprod(A),
  MatMult_C(t(A),A)
)

#######test the Eigen_C function

K[upper_tri(K)] <- K[lower_tri(K)]
temp <- t(A)%*%K%*%A
result <- eigen(temp,symmetric = T)
result2 <- Eigen_C(temp)
all.equal(result$value,result2$value)
microbenchmark(
  eigen(temp,symmetric = T),
  Eigen_C(temp)
)


#####test the colsums function
result <- colSums(temp)
result2 <- colsums(temp)
result3 <- as.vector(ColSum_C(temp))
all.equal(result,result2)
all.equal(result,result3)
microbenchmark(
  colSums(temp),
  colsums(temp),
  as.vector(ColSum_C(temp))
)


##test the Sum_C function
temp <- rnorm(100)
result <- sum(temp)
result2 <- Sum_C(temp)
all.equal(result,result2)
microbenchmark(
  sum(temp),
  Sum_C(temp)
)


##test the MatrixPlus_C function
temp1 <- matrix(rnorm(300*300),300,300)
temp2 <- 0.15
result <- temp2*temp1
result2 <- NumxMatrix_C(temp2,temp1)
all.equal(result,result2)
microbenchmark(
  (temp2*temp1),
  NumxMatrix_C(temp2,temp1)
)

##test the MatrixRowMax_C function
temp1 <- matrix(rnorm(300*300),300,300)


result <- apply(temp1,1,max)
result2 <- MatrixRowMax_C(temp1)
all.equal(result,result2)
microbenchmark(
  apply(temp1,1,max),
  MatrixRowMax_C(temp1)
)



##test the ColSumtwomatrix_C function
temp1 <- matrix(rnorm(300*300),300,300)
temp2 <- matrix(rnorm(300*300),300,300)
result1 <- colSums(temp1)+colSums(temp2)
result2 <- ColSumtwomatrix_C(temp1,temp2)
all.equal(result1,result2)
microbenchmark(
  colSums(temp1)+colSums(temp2),
  ColSumtwomatrix_C(temp1,temp2)
)


##test the ColSumtwomatrix_C function
temp1 <- matrix(rnorm(300*300),300,300)

result1 <- temp1^2
result2 <- Elementwisesquare_C(temp1)
all.equal(result1,result2)
microbenchmark(
  temp1^2,
  Elementwisesquare_C(temp1)
)


##test the crossprod_C function
a <- rnorm(4)
temp <- matrix(rnorm(4*10000),4,10000)
result1 <- (1/a)%*%temp
result2 <- ColSum_C(temp/a)
result3 <- crossprod(1/a,temp)
result4 <- VecMultMat_C(1/a,temp)
all.equal(result1,result2)
all.equal(result1,result3)
all.equal(as.vector(result1),result4)
microbenchmark(
  (1/a)%*%temp ,
  ColSum_C(temp/a),
  crossprod(1/a,temp),
  VecMultMat_C(1/a,temp)
)


###test ifelsetest_C
a <- rnorm(100)
result1 <- ifelse(a<0,0,a)
result2 <- ifelsetest_C(a)
all.equal(result1,result2)
microbenchmark(
  ifelse(a<0,0,a),
  ifelsetest_C(a)
)

### Vecplus_C
a <- rnorm(10000)
b <- rnorm(10000)
result1 <- a+b
result2 <- Vecplus_C(a,b)
all.equal(result1,result2)
microbenchmark(
  a+b,
  Vecplus_C(a,b)
)





Lambdas.test <- Lambdas
mu.test <- mu
w1.test <- w1
w2.test <- w2
nminuspx <- n-px


result1 <- LR0_fixRho_C(Lambdas.test,
                        mu.test ,
             w1.test,
             w2.test,
             nminuspx)

result2 <-  LR0_fixRho
for (i in 1:length.lambda){
  lam = Lambdas.test[i]
  lammu.con <- 1/(1 + lam*mu.test)
  lammu.case <- 1-lammu.con
  Dn = (lammu.con)%*%w1.test+ w2.test
  Nn = (lammu.case)%*%w1.test
  temp = (n-px)*log(1 + Nn/Dn) - Sum_C(log(1 + lam*mu.test))
  result2[,i] = ifelsetest_C(temp)
}

all.equal(result1,result2)
library(microbenchmark)
microbenchmark(LR0_fixRho_C(Lambdas.test,
                            mu.test ,
                            w1.test,
                            w2.test,
                            nminuspx)
,
for (i in 1:length.lambda){
  lam = Lambdas.test[i]
  lammu.con <- 1/(1 + lam*mu.test)
  lammu.case <- 1-lammu.con
  Dn = (lammu.con)%*%w1.test+ w2.test
  Nn = (lammu.case)%*%w1.test
  temp = (n-px)*log(1 + Nn/Dn) - Sum_C(log(1 + lam*mu.test))
  result2[,i] = ifelsetest_C(temp)
})
result1 <- (a*b/(1 + a*b))%*%ww
result2 <- ColSum_C(a*ww*b/(1 + a*b))
all.equal(result1,result2)
microbenchmark(
  (a*b/(1 + a*b))%*%ww,
  ColSum_C(a*ww*b/(1 + a*b))
)
