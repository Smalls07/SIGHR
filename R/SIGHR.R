#' @title A SIGHR Function
#'
#' @description This function allows to leverage side information from other studies to inform about the sparsity structure of the regression coefficients.
#' @param X: the design matrix with confounding variables, SNPs.
#' @param Y: the response variable.
#' @param D: the auxiliary matrix.
#' @param Dtil: list of CMP infor. 1. loc: each row is a location selected from the auxiliary information space. 2. hyper: each row is a set of hyperparameteres for the corresponding location.
#' @param q: the number of confounding (sex, age, BMI... including intercept)
#' @param coef: the intial value of (alpha, beta)
#' @param prob: the intial value inclusion values of beta, of length (p-q)
#' @returns an MCMC sample of all regression coefficients ("coef") including the ones for the intercept and confounders in (1), indicators for all beta's ("indicator.z"), all regression coefficients in the hierarchical logistic regression ("gamma"), and the variance ("sig2") in (1).
#' @export
SIGHR = function(X, Y, D, Dtil, q, coef, prob, iter){
  C = X[,1:q]
  if(q ==1) C = as.matrix(C)
  n = dim(X)[1]
  p = dim(X)[2]  ###number of total covariates with q demographic variables.
  m = dim(D)[2]  ###number of coefficients in the binary model


  ###space holders
  COEF = matrix(-99, nrow = iter, ncol = p)
  Z = matrix(-99, nrow = iter, ncol = p-q)  ###indicator's for if beta should be included or not
  GAMMA = matrix(-99, nrow = iter, ncol = m)	 ###alpha in the simulation code.
  SIG2= rep(-99, iter)

  ###initialization
  d2 = 100 ###prior variance of alpha
  c2 = 1000 ###constant slab variance

  alpha = coef[1:q]
  beta = coef[-c(1:q)]
  prob = prob


  Dstar = Dtil$loc
  apria = Dtil$hyper[1,]
  aprib = Dtil$hyper[2,]


  z = as.numeric(beta != 0)
  omega = rep(1,p-q)
  eta = rep(1,m)
  ni = apria+aprib
  kappa = z - 1/2
  l = apria - (apria+aprib)/2
  sig2 = 1  ###error variance
  gamma = c(-1, 0.1)


  ###log-likelihood with betas integrated out.
  lleval = function(indicator,k){

    if(all(indicator ==0)){
      eval = 0
      return(eval)
    }
    id = which(indicator != 0)
    lid = length(id)
    Xa = X[,id+q]
    sigmainv = t(Xa) %*% Xa / sig2 + diag(lid)/c2
    mu = solve(sigmainv) %*% t(Xa) %*% (Y - C%*%alpha) / sig2
    logdet = determinant(sigmainv/n, T)
    eval = -(1/2) * (logdet[[1]][1]*logdet[[2]] + lid*log(n)) + (1/2) * t(mu) %*% sigmainv %*% mu -(1/2) * log(c2) * indicator[k]
    return(eval)
  }

  for(i in 1:iter){
#    print(i)
    # tic()
    id = sample((q+1):p, (p-q), replace = F)
    for(j in id){
      # print(which(id == j))
      # j = id[1]
      z[j-q] = 1
      f1 = lleval(z,(j-q))
      z[j-q] = 0
      f0 = lleval(z,(j-q))
      f0f1 = exp(f0 - f1)
      probstar = 1/(1+f0f1 * (1-prob[j-q])/prob[j-q])
      # print(c(f0,f1,probstar))
      z[j-q] = rbinom(1,1,as.numeric(probstar))
    }


    ###update beta
    beta[which(z==0)] = 0
    actid = which(z != 0)
    if(length(actid) != 0){
      Xa = X[,(actid+q)]
      sig.beta = solve(t(Xa) %*% Xa/sig2 + diag(length(actid))/c2)
      mu.beta = sig.beta %*% t(Xa) %*% (Y - C%*%alpha)/sig2
      beta[actid] = as.vector(mvrnorm(1, mu.beta, sig.beta))
    }

    ###update kappa
    kappa = z - 1/2

    ###update alpha
    sig.alpha = solve(t(C)%*%C/sig2 + diag(q)/d2)
    mu.alpha = sig.alpha %*% t(C) %*% (Y - X[,-(1:q)] %*% beta)/sig2
    alpha = as.vector(mvrnorm(1, mu.alpha, sig.alpha))

    ###update gamma
    sig.gamma = solve(t(D*omega)%*%D + t(Dstar*eta)%*%Dstar)
    mu.gamma = sig.gamma %*% (t(D) %*% kappa + t(Dstar) %*% l)
    gamma = as.vector(mvrnorm(1,mu.gamma,sig.gamma))


    ###update prob
    prob = 1/(1 + exp(-D %*% gamma))

    ###update omega
    omega = rpg(p-q, 1, D%*%gamma)

    ###update eta
    eta = rpg(m,ni,Dstar%*%gamma)

    ###update sig2
    YCAXB = Y-X%*%c(alpha,beta)
    phi = rgamma(1, shape = n/2 +2, rate = t(YCAXB)%*%YCAXB/2)
    sig2 = 1/phi



    COEF[i,] = c(alpha,beta)
    Z[i,] = z
    GAMMA[i,] = gamma
    SIG2[i] = sig2


  }
  return(list("coef"= COEF, "indicator.z" = Z, "gamma" = GAMMA, "sig2" = SIG2))
}
