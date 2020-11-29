#'@title Estimate L. This only works for z-scores using the ff method.
#'@param Z_hat matrix of z-scores
#'@param R Estimated residual correlation of rows of Z_hat
#'@param fit Flash fit
#'@export
est_L_z <- function(Z_hat, R, fit, opt1 = FALSE, opt2 = TRUE){

  n_var <- nrow(Z_hat)
  n_trait <- ncol(Z_hat)
  n_factor <- ncol(fit$F_hat)
  stopifnot(opt1 | opt2)
  stopifnot(nrow(R) == n_trait & ncol(R) == n_trait)
  stopifnot(all(diag(R) == 1 ))
  stopifnot(matrixcalc::is.positive.definite(R))
  eS <- eigen(R, symmetric=TRUE)

  if(opt1){
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
  }else if(opt2){
    F_total <- fit$fit$loadings.pm[[2]]
    L_total <- fit$fit$loadings.pm[[1]]
    d <- fit$fit$loadings.scale
    ix <- which(colSums(F_total==0) == n_trait)
    F_total <- F_total[,-ix]
    L_total <- L_total[,-ix]
    d <- d[-ix]
    n <- ncol(fit$F_hat)
    F_resid <- F_total[,-(1:n)]
    L_resid <- L_total[, -(1:n)]
    d <- d[-(1:n)]
    X <- L_resid %*% diag(d) %*% t(F_resid)
    Sigma <- diag(fit$fit$residuals.sd^2) + (t(X) %*% X)/(nrow(L_resid))
  }
  L_est_var <-  H%*% Sigma %*% t(H)
  L_est_se <- t( sqrt(diag(L_est_var)) * t(matrix(1, nrow=n_var, ncol=n_factor)))

  return(list(L_est = L_est, L_est_se = L_est_se))
}
