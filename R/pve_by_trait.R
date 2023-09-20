
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
