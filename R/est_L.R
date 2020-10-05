#'@export
est_L <- function(B_hat, S_hat, N, tau, adjust=TRUE, fit, ignore_tau =FALSE){

  n_var <- nrow(B_hat)
  n_trait <- ncol(B_hat)
  n_factor <- ncol(fit$F_hat)
  if(!ignore_tau){
    if(missing(tau)){stop("Please provide tau.")}
    stopifnot(length(tau)==n_trait)
  }else{
    tau <- rep(1, n_trait)
  }
  if(adjust){
    B_tilde = t( (1/sqrt(N)) *t(B_hat/S_hat))
    S2_tilde = t( (1/N) * t(matrix(1, nrow=n_var, ncol=n_trait)))
  }else{
    B_tilde = B_hat
    S2_tilde = S_hat^2
  }

  H <- with(fit, solve(t(F_hat) %*% F_hat) %*% t(F_hat))
  L_est <- H %*% t(B_hat) %>% t()
  L_est_se <- map(seq(n_var), function(i){
    V <- H%*% diag(S2_tilde[i,] + 1/tau) %*% t(H)
    sqrt(diag(V))
  }) %>%
  do.call(rbind, .)

  return(list(L_est = L_est, L_est_se = L_est_se))
}
