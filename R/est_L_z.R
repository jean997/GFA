#'@title Estimate L. This only works for z-scores using the ff method.
#'@param Z_hat matrix of z-scores
#'@param R Estimated residual correlation of rows of Z_hat
#'@param fit Flash fit
#'@export
est_L_z <- function(Z_hat, R, fit){

  n_var <- nrow(Z_hat)
  n_trait <- ncol(Z_hat)
  n_factor <- ncol(fit$F_hat)

  stopifnot(nrow(R) == n_trait & ncol(R) == n_trait)
  stopifnot(all(diag(R) == 1 ))
  stopifnot(matrixcalc::is.positive.definite(R))
  eS <- eigen(R, symmetric=TRUE)

  if(all(R == diag(1, n_trait))){
    R_is_identity <- TRUE
  }else{
    R_is_identity <- FALSE
    lambda_min <- min(eS$values)
    V <- eS$vectors[, -n_trait]
    W <- V %*% diag(sqrt(eS$values[-n_trait]-lambda_min)) %*% t(V)
  }

  H <- with(fit, solve(t(F_hat) %*% F_hat) %*% t(F_hat))
  L_est <- H %*% t(Z_hat) %>% t()
  Sigma <- diag(fit$fit$residuals.sd^2) + W
  L_est_var <-  H%*% Sigma %*% t(H)
  L_est_se <- t( sqrt(diag(L_est_var)) * t(matrix(1, nrow=n_var, ncol=n_factor)))

  return(list(L_est = L_est, L_est_se = L_est_se))
}
