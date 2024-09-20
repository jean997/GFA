loadings_gls <- function(X, S, R, F_hat){
  ntrait <- ncol(X)
  nvar <- nrow(X)
  nfactor <- ncol(F_hat)
  stopifnot(nrow(F_hat) == ntrait)
  stopifnot(nrow(S) == nvar & ncol(S) == ntrait)

  if(missing(R)){
    R <- diag(ntrait, nrow = ntrait)
  }
  stopifnot(nrow(R) == ntrait & ncol(R) == ntrait)
  stopifnot(isSymmetric(R))
  stopifnot(all(eigen(R)$values > 0))

  Rinv <- solve(R)
  if(all(S == 1)){
    ftf_inv <- solve(t(F_hat) %*% Rinv %*% F_hat)
    ftx <- t(F_hat) %*% Rinv %*% t(X)
    L <- t(ftf_inv %*% ftx)
    S <- t(matrix(1, nrow = nfactor, ncol = nvar) * sqrt(diag(ftf_inv)))
  }else{
    LS <- sapply(1:nvar, function(i){
      s <- S[i,]
      ftf_inv <- solve(t(F_hat) %*% diag(1/s) %*% Rinv %*% diag(1/s) %*% F_hat)
      ftx <- t(F_hat) %*% diag(1/s) %*% Rinv %*% diag(1/s) %*% t(X[i,,drop=F])
      c(ftf_inv %*% ftx, sqrt(diag(ftf_inv)))
    }) %>% matrix(nrow = nvar, byrow = TRUE)
    L <- LS[, 1:nfactor]
    S <- LS[, -(1:nfactor)]
  }
  return(list("L" = L, "S" = S))
}


#'@export
gfa_loadings_gls <- function(beta_hat, S, fit){

  if(!fit$mode == "z-score"){
    stop("mode must be z-score to use this function.\n")
  }

  X <- beta_hat/S
  if(is.null(fit$R)){
    myR <- diag(fit$fit$residual_sd^2)
  }else if(fit$method == "fixed_factors"){
    myR <- fit$R + diag(fit$fit$residuals_sd^2)
  }else if(fit$method == "noR"){
    myR <- fit$R + diag(fit$fit$residuals_sd^2-1)
  }
  myS <- matrix(1, nrow = nrow(X), ncol = ncol(X))
  myF <- fit$F_hat*fit$scale ## put scale back in because we are using z-scores
  ret <- loadings_gls(X = X, S = myS, R = myR, F_hat = myF)
  ret$P <- 2*pnorm(-abs(ret$L/ret$S))
  return(ret)
}
