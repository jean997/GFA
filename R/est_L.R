#'@title Estimate L
#'@param B_hat matrix of effect estimates
#'@param S_hat matrix of standard errors
#'@param R Estimated residual correlation of rows of B_hat
#'@param N Sample size
#'@param tau vector length n traits of per-trait precision (1/variance) of theta
#'@param fit Flash fit
#'@param adjust. Whether to adjust estimates to use standardized effects.
#'@export
est_L <- function(B_hat, S_hat, R, N, tau, fit,  adjust=TRUE, tol=1e-15){

  n_var <- nrow(B_hat)
  n_trait <- ncol(B_hat)
  n_factor <- ncol(fit$F_hat)

  stopifnot(length(tau)==n_trait)
  stopifnot(nrow(R) == n_trait & ncol(R) == n_trait)
  stopifnot(all(diag(R) == 1 ))
  stopifnot(matrixcalc::is.positive.definite(R))


  if(all(R == diag(1, n_trait))){
    R_is_identity <- TRUE
  }else{
    R_is_identity <- FALSE
  }

  if(adjust){
    B_tilde = t( (1/sqrt(N)) *t(B_hat/S_hat))
    S_tilde = t( (1/sqrt(N)) * t(matrix(1, nrow=n_var, ncol=n_trait)))
    all_snps_same_S <- TRUE
  }else{
    B_tilde = B_hat
    S_tilde = S_hat
    if(all(apply(S_tilde, 2, function(x){all(x==x[1])}))){
      all_snps_same_S <- TRUE
    }else{
      all_snps_same_S <- FALSE
    }
  }

  H <- with(fit, solve(t(F_hat) %*% F_hat) %*% t(F_hat))
  L_est <- H %*% t(B_hat) %>% t()


  if(all_snps_same_S){
    S <- diag(S_tilde[1,] + 1/sqrt(tau))
    V <- H%*% S %*% R %*% S %*% t(H)
    L_est_se <- t( sqrt(diag(V)), t(matrix(1, nrow=n_var, ncol=n_trait)))
  }else if(R_is_identity){
    L_est_se <- map(seq(n_var), function(i){
                    V <- H%*% diag(S_tilde[i,]^2 + 1/tau) %*% t(H)
                    sqrt(diag(V))
                }) %>%
                do.call(rbind, .)
  }else{
    L_est_se <- map(seq(n_var), function(i){
                    S <- diag(S_tilde[i,] + 1/sqrt(tau))
                    V <- H%*% S %*% R %*% S %*% t(H)
                    sqrt(diag(V))
                    }) %>%
                do.call(rbind, .)
  }

  return(list(L_est = L_est, L_est_se = L_est_se))
}
