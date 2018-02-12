
#'
#' AKA_C <- function(A,K){
#'   AV <- as.vector(A)
#'   KV <- as.vector(K)
#'   A.row <- nrow(A)
#'   A.col <- ncol(A)
#'   outV <- rep(0,A.col^2)
#'   #dyn.load("tAKA.so")
#'   temp <- .C("tAKA",AV,KV,A.row,A.col,outV=outV)
#'   result <- matrix(temp[[5]],A.col,A.col)
#'   return(result)
#' }
#' #t(A)%*%K%*%A
