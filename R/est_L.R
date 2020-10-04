#'@export
est_L <- function(B_hat, S_hat, N, adjust=TRUE, fit, ignore_tau =TRUE){
  stopifnot(ignore_tau)
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

  H <- with(fit, solve(t(F_hat) %*% F_hat) %*% t(F_hat))
  L_est <- H %*% t(B_hat) %>% t()
  if(ignore_tau){
    L_est_var <- map(seq(n_var), function(i){
                          V <- H%*% diag(S_hat[i,]) %*% t(H)
                          sqrt(diag(V))
                        }) %>%
                  do.call(rbind, .)
  }else{
    tau <- flashier:::get.tau(fit$fit$flash.fit)
    T2_tilde = t( (1/tau) * t(matrix(1, nrow=n_var, ncol=n_trait)))
  }
  return(list(L_est = L_est, L_est_se = L_est_se))
}
