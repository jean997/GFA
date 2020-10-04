#'@export
est_L <- function(B_hat, S_hat, N, adjust=TRUE, fit){
  n_var <- nrow(B_hat)
  n_trait <- ncol(B_hat)
  n_factor <- ncol(fit$F_hat)
  if(adjust){
    B_tilde = t( (1/sqrt(N)) *t(B_hat/S_hat))
    S2_tilde = t( (1/N) * t(matrix(1, nrow=n_var, ncol=n_trait)))
  }else{
    B_tilde = B_hat
    S2_tilde = S_hat^2
  }
  tau <- flashier:::get.tau(fit$fit$flash.fit)
  T2_tilde = t( (1/tau) * t(matrix(1, nrow=n_var, ncol=n_trait)))

  H <- solve(t(F_hat) %*% F_hat) %*% t(F_hat)
  L_est <- H %*% t(B_hat) %>% t()
  S <- diag(S2_tilde[1,] + T2_tilde[1,])
  L_est_var <- H%*% S %*% t(H)
  L_est_se <- t(sqrt(diag(L_est_var))*t(matrix(1, n_var,n_factor)))
  return(list(L_est = L_est, L_est_se = L_est_se))
}
