
#' Title
#'
#' @param A
#' @param K
#'
#' @return
#' @export
#'
#' @examples
tAKA <- function(A,K){
  AV <- as.vector(A)
  KV <- as.vector(K)
  A.row <- nrow(A)
  A.col <- ncol(A)
  out <- rep(0,A.col^2)
  temp <- .C("TAKA",AV,KV,A.row,A.col,out=out)
  result <- as.matrix(temp[[5]],A.col,A.col)
}
