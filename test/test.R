
y=y; X = cbind(X, E);K1 = G %*% t(G); K2 = (G*E) %*% t(G * E) ; N = 10000; length.rho = 200; length.lambda = 21

X = X[,c(1,2)]

  method = "ML"
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

  eP   = eigen(P0, symmetric = T)
  A    = eP$vector[,eP$values > 1e-10]
  # invP = ginv(P0)
  # A2   = eP$vectors %*% diag(eP$values)
  eK1  = eigen(K1, symm = T)
  wK1  = which(eK1$values  > 1e-10)
  # phi1 = t(t(eK1$vectors[,wK1])*sqrt(eK1$values[wK1]))
  if (length(wK1) == 1){
    phi1 = eK1$vectors[,wK1] * sqrt(eK1$values[wK1])
  }else{
    phi1 = t(t(eK1$vectors[,wK1])*sqrt(eK1$values[wK1]))
  }

  eK2  = eigen(K2, symm = T)
  wK2  = which(eK2$values  > 1e-10)

  if (length(wK2) == 1){
    phi2 = eK2$vectors[,wK2] * sqrt(eK2$values[wK2])
  }else{

    phi2 = t(t(eK2$vectors[,wK2])*sqrt(eK2$values[wK2]))
  }

  group= rep(1,n)
  fit1 = lme(y~X, random = list(group=pdIdent(~-1+phi1), group = pdIdent(~-1+phi2)), method = "ML") #  Default = REML
  fit0 = lm(y~X)
  LR = max(0, 2*(logLik(fit1, REML = F) -logLik(fit0, REML = F)))

  # if (LR <= 0){
  #   p.dir = 1;
  #   p.au1=1; p.aud = 1
  #
  # }else{
  #   #For the first kernel
    LR0_allRho = matrix(NA, N, length.rho)
    set.seed(123)
    w = matrix(rnorm(N*(n-px)), n-px,N)
    LR0_fixRho = matrix(NA, N, length.lambda)
    # The baseline W1
    rho = 0   # K = rho*K1 + (1-rho)*K2 = K2
    K = K2
    k = length(wK2)
    xi = eK2$values[wK2]

    mu = eigen(t(phi2) %*% P0%*% phi2, symmetric = T,  only.values = T)$values
    # mu = eigen(P0 %*% K2)$values

    mu = mu/max(mu, xi)
    xi = xi/max(mu,xi)
    AKA = t(A) %*% K %*% A
    eV1 = eigen(AKA, symmetric = T)
    # eV1= eV
    # U1 = eV1$vectors
    w1 = (w^2)[1:k,]
    w2 = colSums((w^2)[-(1:k),])

    if (length(mu) < k){mu = c(mu,rep(0, k - length(mu)))}
    if (length(xi) < k){xi = c(xi,rep(0, k - length(xi)))}

    for (i in 1:length.lambda){
      lam = Lambdas[i]
      Dn = colSums(w1/(1 + lam*mu))+ w2
      Nn = colSums(lam*w1*mu/(1 + lam*mu))
      temp = (n-px)*log(1 + Nn/Dn) - sum(log(1 + lam*xi))
      LR0_fixRho[,i] = ifelse(temp < 0, 0, temp)
    }
    LR0_allRho[,1] = apply(LR0_fixRho, 1, max)

   # for (j in 2:length.rho){
    j <- 2
      rho = all_rho[j]
      LR0_fixRho = matrix(NA, N, length.lambda)
      K = rho*K1 +(1-rho)*K2
      eK = eigen(K, symmetric = T)
      wK = which(eK$values > 1e-10)
      k = length(wK)
      xi = eK$values[wK]

      phi= t(t(eK$vectors[,wK]) * sqrt(eK$values[wK]) )
      mu = eigen(t(phi) %*% P0%*% phi, symmetric = T,  only.values = T)$values

      mu = mu/max(mu, xi)
      xi = xi/max(mu,xi)


      AKA= t(A) %*% K %*% A
      eV = eigen(AKA, symmetric = T)
      u2 = eV$vector
      # W2 = A %*% eV$vectors

      uu = t(eV$vector)%*% eV1$vector
      ww = uu %*% w
      ww.double.old = ww^2
      w1 = (ww^2)[1:k,]
      w2 = colSums((ww^2)[-(1:k),])
      if (length(mu) < k){mu = c(mu,rep(0, k - length(mu)))}
      if (length(xi) < k){xi = c(xi,rep(0, k - length(xi)))}


      for (i in 1:length.lambda){
        lam = Lambdas[i]
        Dn = colSums(w1/(1 + lam*mu))+ w2
        Nn = colSums(lam*w1*mu/(1 + lam*mu))
        temp = (n-px)*log(1 + Nn/Dn) - sum(log(1 + lam*xi))
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
  out = list(p.dir = p.dir,p.au1=p.au1,p.aud = p.aud, LR = LR)



