#'@title GFA generalized least squares loadings
#'@description Compute GLS estimates and p-values for loadings
#'@param beta_hat Variants by traits matrix of association estimates. Traits
#'should be in the same order as used to produce fitted gfa object. The variant set may be different than
#'the variant set used to produce the fitted gfa object.
#'@param S Variants by traits matrix of standard errors. Variants and traits should match variants and traits in beta_hat.
#'@param fit Object produced by gfa_ft().
#'@return A list with elements L (estimated loadings), S (standard errors for loadings), and P (p-values). All variant by factor matrices.
#'@export
gfa_loadings_gls <- function(beta_hat, S, fit){

  #if(!fit$mode == "z-score"){
  #  stop("mode must be z-score to use this function.\n")
  #}

  X <- beta_hat/S
  if(is.null(fit$R) | fit$method == "noR"){
    myR <- diag(fit$fit$residuals_sd^2)
  }else if(fit$method == "fixed_factors"){
    myR <- fit$R + diag(fit$fit$residuals_sd^2)
  }

  myS <- matrix(1, nrow = nrow(X), ncol = ncol(X))
  myF <- fit$F_hat*fit$scale ## put scale back in because we are using z-scores
  ret <- loadings_gls(X = X, S = myS, R = myR, F_hat = myF)
  ret$P <- 2*pnorm(-abs(ret$L/ret$S))
  return(ret)
}





loadings_gls <- function(X, S, R, F_hat){
  ntrait <- ncol(X)
  nvar <- nrow(X)
  nfactor <- ncol(F_hat)
  stopifnot(nrow(F_hat) == ntrait)
  stopifnot(nrow(S) == nvar & ncol(S) == ntrait)

  if(missing(R)){
    R <- diag(1, nrow = ntrait)
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

