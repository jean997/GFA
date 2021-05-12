
# Should return a matrix that is traits by factors
pve_by_trait <- function(Lhat, Fhat, tau){
  nvar <- nrow(Lhat)
  ntrait <- nrow(Fhat)
  k <- ncol(Lhat)
  stopifnot(ncol(Fhat) == k)
  if(class(tau)=="numeric") tau <- matrix(tau, nrow=nvar, ncol=ntrait)
  stopifnot(nrow(tau) == nvar & ncol(tau) == ntrait)
  sj <- purrr::map_dfc(seq(k), function(kk){
    B <- Lhat[,kk] %*% t(Fhat[,kk])
    colSums(B^2)
  })
  tau_sums <- colSums(1/tau)
  s_sums <- rowSums(sj)
  sums <- matrix(rep(s_sums + tau_sums, k), byrow = F)
  pve_j <- sj/sums
}

#'@export
pve_by_factor2 <- function(fit, S_hat){
  stopifnot(class(fit$residuals.sd) == "numeric")
  nvar <- nrow(fit$loadings.pm[[1]])
  ntrait <- nrow(fit$loadings.pm[[2]])
  k <- ncol(fit$loadings.pm[[1]])
  if(missing(S_hat)){
    S_hat = matrix(1, nrow = nvar, ncol = ntrait)
  }
  stopifnot(nrow(S_hat) == nvar & ncol(S_hat) == ntrait)
  sj <- purrr::map_dfc(seq(k), function(kk){
    B <- fit$loadings.pm[[1]][,kk] %*% t(fit$loadings.pm[[2]][,kk])
    B <- B*(fit$loadings.scale[kk]^2)*(S_hat^2)
    colSums(B^2)
  })
  trait_sums <- rowSums(sj)
  theta_var <- (S_hat^2)*matrix(rep(fit$residuals.sd^2, nvar), byrow=TRUE, nrow = nvar)
  theta_trait_sums <- colSums(theta_var)
  total_var <- trait_sums + theta_trait_sums
  pve_by_trait <- sj/total_var
  ix <- which(fit$loadings.scale ==0)
  pve_by_trait <- pve_by_trait[,-ix]
  return(pve_by_trait)
}

