library(Utility)
A.row <- 300
A.col <- 297
A <- matrix(rnorm(A.row*A.col),A.row,A.col)
K <- matrix(rnorm(A.row*A.row),A.row,A.row)
result <- t(A)%*%K%*%A
result2 <- tAKA(A,K)
