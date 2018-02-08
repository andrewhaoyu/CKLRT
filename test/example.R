library(devtools)
install_github("andrewhaoyu/Utility")
try(library(mgcv)) ; try(library(MASS)); library(nlme);
library(compiler);library(Rcpp);library(RcppEigen)
library(Utility)
n = 1000
X = rnorm(n)

p = 2
G = runif(n*p)< 0.5
G = G + runif(n*p) < 0.5
G = matrix(G, n,p)
E = (runif(n) < 0.5)^2
y = rnorm(n) + G[,1] * 0.3

time1 <- proc.time()
omniRLRT(y, X = cbind(X, E),K1 = G %*% t(G), K2 = (G*E) %*% t(G * E) , N = 10000, length.rho = 200, length.lambda = 21)
time1 <- proc.time()-time1
time2 <- proc.time()
omniRLRT_fast(y, X = cbind(X, E),K1 = G %*% t(G), K2 = (G*E) %*% t(G * E) , N = 10000, length.rho = 200, length.lambda = 21)
time2 <- proc.time()-time2
print(time1)
print(time2)
# microbenchmark(
#   omniRLRT(y, X = cbind(X, E),K1 = G %*% t(G), K2 = (G*E) %*% t(G * E) , N = 10000, length.rho = 200, length.lambda = 21),
#   omniRLRT_fast(y, X = cbind(X, E),K1 = G %*% t(G), K2 = (G*E) %*% t(G * E) , N = 10000, length.rho = 200, length.lambda = 21)
# )
