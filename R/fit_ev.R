
#'@title Fit with eigenvector transformation
#'@param B_hat Matrix of (non-standardized) effect estimates SNPs by traits
#'@param S_hat Matrix of standard errors SNPs by traits
#'@param N Vector of sample sizes length equal to number of traits
#'@param R Estimated residual correlation matrix
#'@param kmax Maximum number of factors
#'@param zero_thresh Threshold for setting eigenvalues of R to zero
#'@return A list with elements fit, B_hat, L_hat, F_hat
#'@export
fit_ev <- function(B_hat, S_hat, N, R, kmax=100, zero_thresh = 1e-15){

  n_var <- nrow(B_hat)
  n_trait <- ncol(B_hat)
  stopifnot(nrow(S_hat) == n_var & ncol(S_hat) == n_trait)
  stopifnot(length(N) == n_trait)

  R_eig <- eigen(R)
  R_eig$values[abs(R_eig$values) < zero_thresh] <- 0
  if(any(R_eig$values < 0))stop("R is not psd")
  d <- R_eig$values
  V <- R_eig$vectors

  B_tilde = t( (1/sqrt(N)) *t(B_hat/S_hat))
  S_tilde = t( (1/sqrt(N)) * t(matrix(1, nrow=n_var, ncol=n_trait)))

  B_tilde_tilde <- B_tilde %*% V
  S_tilde_tilde <- t( (sqrt(d)/sqrt(N)) * t(matrix(1, nrow=n_var, ncol=n_trait)))

  fit <- flash.init(data=B_tilde_tilde, S = S_tilde_tilde,  var.type=2) %>%
    flash.add.greedy(Kmax = kmax, init.fn = init.fn.softImpute) %>%
    flash.backfit() %>%
    flash.nullcheck();

  if(fit$n.factors > 0){
    F_hat <- V %*% fit$loadings.pm[[2]]
    L_hat <- fit$loadings.pm[[1]]
    B_hat <- fitted(fit) %*% t(V)
  }else{
    B_hat <- F_hat  <- L_hat <- NULL;
  }
  ret <- list(fit=fit, B_hat = B_hat, L_hat = L_hat, F_hat = F_hat)
  return(ret)
}
