#'@title Estimate L. This only works for z-scores using the ff method.
#'@param Z_hat matrix of z-scores
#'@param R Estimated residual correlation of rows of Z_hat
#'@param fit Flash fit
#'@export
est_L_z <- function(Z_hat, R, F_hat){

  n_var <- nrow(Z_hat)
  n_trait <- ncol(Z_hat)
  n_factor <- ncol(F_hat)

  stopifnot(nrow(R) == n_trait & ncol(R) == n_trait)
  #stopifnot(all(diag(R) == 1 ))
  stopifnot(matrixcalc::is.positive.definite(R, tol= 1e-10))
  Rinv <- solve(R)

  A <- solve(t(F_hat) %*% Rinv %*% F_hat)
  H <- A %*% (t(F_hat) %*% Rinv)
  L_est <- H %*% t(Z_hat) %>% t()
  #L_est_var <-  H%*% Sigma %*% t(H)

  L_est_se <- t( sqrt(diag(A)) * t(matrix(1, nrow=n_var, ncol=n_factor)))

  return(list(L_est = L_est, L_est_se = L_est_se))
}
