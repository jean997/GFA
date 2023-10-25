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
