#' Title
#'
#' @param y
#' @param X
#' @param K1
#' @param K2
#' @param N
#' @param length.rho
#' @param length.lambda
#'
#' @return
#' @export
#'
#' @examples
omniRLRT_fast = function(y, X,K1, K2, N = 10000, length.rho = 200, length.lambda = 21){
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
  XX   = MatMult_C(t(X1),X1)
  P0   = diag(n)- MatMult_C(MatMult_C(X1,ginv(XX)),t(X1))

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

  if (LR <= 0){
    p.dir =p.au1= p.aud = 1
  }else{

    #For the first kernel
    LR0_allRho = matrix(NA, N, length.rho)
    w = matrix(rnorm(N*(n-px)), n-px,N)
    LR0_fixRho = matrix(NA, N, length.lambda)

    rho = 0
    K = K2
    k = length(wK2)
    xi = eK2$values[wK2]
    AKA = MatMult_C(MatMult_C(t(A),K),A)
    eV = Eigen_C(AKA)
    #eV = Eigen_C(AKA)
    U_1 = eV$vectors
    mu = eV$values[eV$values > 1e-10]


    mu = mu/max(mu, xi)
    xi = xi/max(mu,xi)

    # original U_1 in this case
    W1 = MatMult_C(A,U_1)
    w.double <- w^2
    w1 = w.double[1:k,]
    w2 = ColSum_C((w.double)[-(1:k),])

    if (length(mu) < k){mu = c(mu,rep(0, k - length(mu)))}
    if (length(xi) < k){xi = c(xi,rep(0, k - length(xi)))}


    LR0_fixRho <- LR0_fixRho_C(Lambdas,
                 mu,
                 w1,
                 w2,
                 n-px)
    # for (i in 1:length.lambda){
    #   lam = Lambdas[i]
    #   Dn = (1/(1 + lam*mu))%*%w1+ w2
    #   Nn = (lam*mu/(1 + lam*mu))%*%w1
    #   temp = (n-px)*log(1 + Nn/Dn) - Sum_C(log(1 + lam*mu))
    #   LR0_fixRho[,i] = ifelse(temp < 0, 0, temp)
    # }
     LR0_allRho[,1] = MatrixRowMax_C(LR0_fixRho)

     LR0_allRho <- doubleloop(K1,
                                   K2,
                                   P0,
                                   A,
                                   U_1,
                                   w,
                                   Lambdas,
                                   n-px,
                                   all_rho,
                                   LR0_allRho)
    LR0 = MatrixRowMax_C(LR0_allRho)
    LR0 = ifelse(LR0 > 0, LR0, 0)
    p.dir = mean(LR < LR0)
    p.au1 = getp_au1(null = LR0, LR = LR)$p
    p.aud= getp_aud_estimate_pi_first(null = LR0, LR = LR)$p
  }
  out = list(p.dir = p.dir,p.au1 = p.au1,p.aud = p.aud, LR = LR)
  return(out)
}


