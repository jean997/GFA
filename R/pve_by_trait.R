
# Should return a matrix that is traits by factors
pve_by_trait <- function(Lhat, Fhat, tau, sample_size){
  nvar <- nrow(Lhat)
  ntrait <- nrow(Fhat)
  k <- ncol(Lhat)
  stopifnot(ncol(Fhat) == k)
  if(missing(sample_size)) sample_size <- rep(1, ntrait)
  stopifnot(length(sample_size) == ntrait)
  if(class(tau)=="numeric") tau <- matrix(tau, nrow=nvar, ncol=ntrait, byrow = TRUE)
  stopifnot(nrow(tau) == nvar & ncol(tau) == ntrait)
  # variance explained per factor
  sj <- purrr::map_dfc(seq(k), function(kk){
    B <- Lhat[,kk] %*% t(Fhat[,kk]) ## total effects on original (usually z-score scale)
    B_e <- t(t(B)/sqrt(sample_size)) ## standardized effect scale effects
    colSums(B_e^2)
  })
  tau <- t(t(tau)*sample_size)
  tau_sums <- colSums(1/tau) # per-trait variance explained by error
  s_sums <- rowSums(sj) # per-trait variance explained by factors
  sums <- s_sums + tau_sums
  pve_j <- sj/sums
}

