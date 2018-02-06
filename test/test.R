
y=y; X = cbind(X, E);K1 = G %*% t(G); K2 = (G*E) %*% t(G * E) ; N = 10000; length.rho = 200; length.lambda = 21





method = "REML"
Lambdas = exp(seq(from = -12, to = 12, length.out = length.lambda))
all_rho = seq(from = 0,to = 1, length.out = length.rho)
n  = length(y)
if (is.null(X)){
  X1 = matrix(1, nrow=n)
  px = 1
}else{
  X1 = cbind(1,X)
  px = ncol(X1)
}
XX   = t(X1) %*% X1
P0   = diag(n)- X1 %*% ginv(XX) %*% t(X1)

eP   = Eigen_C(P0)
A    = eP$vector[,eP$values > 1e-10]
# invP = ginv(P0)
# A2   = eP$vectors %*% diag(eP$values)
eK1  = Eigen_C(K1)
wK1  = which(eK1$values  > 1e-10)
# phi1 = t(t(eK1$vectors[,wK1])*sqrt(eK1$values[wK1]))
if (length(wK1) == 1){
  phi1 = eK1$vectors[,wK1] * sqrt(eK1$values[wK1])
}else{
  phi1 = t(t(eK1$vectors[,wK1])*sqrt(eK1$values[wK1])) # this actually works for all of them
}

eK2  = Eigen_C(K2)
wK2  = which(eK2$values  > 1e-10)

if (length(wK2) == 1){
  phi2 = eK2$vectors[,wK2] * sqrt(eK2$values[wK2])
}else{
  phi2 = t(t(eK2$vectors[,wK2])*sqrt(eK2$values[wK2]))
}
# if wK1 and wK2 are 1, what would this happen to others?

group= rep(1,n)
fit1 = lme(y~X, random = list(group=pdIdent(~-1+phi1), group = pdIdent(~-1+phi2))) #  Default = REML
fit0 = lm(y~X)
LR = max(0, 2*(logLik(fit1, REML = T) -logLik(fit0, REML = T)))

# if (LR <= 0){
#   p.dir =p.au1= p.aud = 1
# }else{

  #For the first kernel
  LR0_allRho = matrix(NA, N, length.rho)
  w = matrix(rnorm(N*(n-px)), n-px,N)
  LR0_fixRho = matrix(NA, N, length.lambda)

  rho = 0
  K = K2
  k = length(wK2)
  xi = eK2$values[wK2]
  xi.2 = xi


  AKA = t(A) %*% K %*% A
  eV = Eigen_C(AKA)
  U_1 = eV$vectors
  mu = eV$values[eV$values > 1e-10]
  eV.2 = eigen(AKA)
  U_1.2 = eV.2$vectors
  mu.2 = eV.2$values[eV.2$values > 1e-10]
  # AKA = t(A) %*% K %*% A
  # eV = eigen(AKA)
  # U_1 = eV$vectors
  # mu = eV$values[eV$values > 1e-10]

  mu = mu/max(mu, xi)
  xi = xi/max(mu,xi)
  mu.2 = mu.2/max(mu.2, xi.2)
  xi.2 = xi.2/max(mu,xi.2)
  # original U_1 in this case
  W1 = A %*% U_1
  W1.2 = length(A%*%U_1.2)
  w1 = (w^2)[1:k,]
  w2 = colSums((w^2)[-(1:k),])

  if (length(mu) < k){mu = c(mu,rep(0, k - length(mu)))}
  if (length(xi) < k){xi = c(xi,rep(0, k - length(xi)))}

  for (i in 1:length.lambda){
    lam = Lambdas[i]
    Dn = colSums(w1/(1 + lam*mu))+ w2
    Nn = colSums(lam*w1*mu/(1 + lam*mu))
    temp = (n-px)*log(1 + Nn/Dn) - sum(log(1 + lam*mu))
    LR0_fixRho[,i] = ifelse(temp < 0, 0, temp)
  }
  LR0_allRho[,1] = apply(LR0_fixRho, 1, max)

  for (j in 2:length.rho){
    rho= all_rho[j]
    LR0_fixRho = matrix(NA, N, length.lambda)
    K  = rho*K1 +(1-rho)*K2
    eK = eigen(K, symmetric = T)
    wK = which(eK$values > 1e-10)
    k  = length(wK)
    xi = eK$values[wK]
    # phi= eK$vectors[,wK] %*% diag(sqrt(eK$values[wK]))
    phi= t(t(eK$vectors[,wK]) * sqrt(eK$values[wK]) )
    mu = eigen(t(phi) %*% P0%*% phi, symmetric = T,  only.values = T)$values

    # mu should be the same as eigen AKA?
    mu = mu/max(mu, xi)
    xi = xi/max(mu,xi)
    # A is eigen vector of P_0
    AKA = t(A) %*% K %*% A
    eV = eigen(AKA, symmetric = T)
    W2 = A %*% eV$vectors
    U2 = eV$vectors
    ww = t(U2) %*% U_1 %*% w

    # ww = t(W2) %*% invP %*% W1  %*% w  # ww is the rotated version of w.
    # By the proof, it is also ww = t(U_2) %*% U_1 %*% w
    w1 = (ww^2)[1:k,]
    w2 = colSums((ww^2)[-(1:k),])

    if (length(mu) < k){mu = c(mu,rep(0, k - length(mu)))}
    if (length(xi) < k){xi = c(xi,rep(0, k - length(xi)))}

    for (i in 1:length.lambda){
      lam = Lambdas[i]
      Dn = colSums(w1/(1 + lam*mu))+ w2
      Nn = colSums(lam*w1*mu/(1 + lam*mu))
      temp = (n-px)*log(1 + Nn/Dn) - sum(log(1 + lam*mu))
      LR0_fixRho[,i] = ifelse(temp < 0, 0, temp)
    }
    LR0_allRho[,j] = apply(LR0_fixRho, 1, max)
  }
  LR0 = apply(LR0_allRho, 1, max)
  LR0 = ifelse(LR0 > 0, LR0, 0)
  p.dir = mean(LR < LR0)
  p.au1 = getp_au1(null = LR0, LR = LR)$p
  p.aud= getp_aud_estimate_pi_first(null = LR0, LR = LR)$p
}
out = list(p.dir = p.dir,p.au1 = p.au1,p.aud = p.aud, LR = LR)














library(devtools)
#install_github("andrewhaoyu/Utility")
library(Utility)
A.row <- 300
A.col <- 297
A <- matrix(rnorm(A.row*A.col),A.row,A.col)
K <- matrix(rnorm(A.row*A.row),A.row,A.row)

#########t(A)%*%K%*%A faster
result <- t(A)%*%K%*%A

result <- matrix(rnorm(300*300),300,300)
result[lower.tri(result)] <- result[upper.tri(result)]
det(result)

try <- eigen(result,symmetric=T)
try.2 <- Eigen_C(result)
try
try.2
try.2$values
P0 <- result

eP   = eigen(P0, symmetric = T)
A    = eP$vector[,eP$values > 1e-10]
# invP = ginv(P0)
# A2   = eP$vectors %*% diag(eP$values)
eK1  = eigen(K1, symm = T)
wK1  = which(eK1$values  > 1e-10)









try2 <- Eigen_C_value(result)
library(microbenchmark)
eigen(result,symmetric = T)
Eigen_C(result)
microbenchmark(
  eigen(result,symmetric = T),
  Eigen_C(result)
)

result2 <- AKA_C(A,K)

#########


