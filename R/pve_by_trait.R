
# Should return a matrix that is traits by factors
pve_by_trait <- function(fit, sample_size){
  nvar <- nrow(fit$L_hat)
  ntrait <- nrow(fit$F_hat)

  if(missing(sample_size)) sample_size <- rep(1, ntrait)
  stopifnot(length(sample_size) == ntrait)

  ef2 <- flashier:::get.EF2(fit$fit$flash_fit)
  nf <- ncol(ef2[[1]])
  est_ix <- setdiff(seq(nf), fit$fixed_ix) |> sort()

  # variance explained per trait/factor
  sj <- sapply(est_ix, function(kk){
    Vk <- ef2[[1]][,kk] %*% t(ef2[[2]][,kk]) ## total effects on original (usually z-score scale)
     ## standardized effect scale total variance
    colSums(Vk)/sample_size
  }) |> matrix(nrow = ntrait, byrow = FALSE)

  # variance from fixed factors (this is part of error)
  fj <- sapply(fit$fixed_ix, function(kk){
    Vk <- ef2[[1]][,kk] %*% t(ef2[[2]][,kk]) ## total effects on original (usually z-score scale)
    ## standardized effect scale total variance
    colSums(Vk)/sample_size
  })|> matrix(nrow = ntrait, byrow = FALSE)

  tau <- flashier:::get.tau(fit$fit$flash_fit)
  stopifnot(length(tau) == ntrait)
  var_from_tau <- nvar/tau
  var_from_tau <- var_from_tau/sample_size

  tot_by_trait <- rowSums(sj) + rowSums(fj) + var_from_tau

  # error component of tau
  fixed_tau <- flash_fit_get_fixed_tau(fit$fit$flash_fit)

  error_var <- rowSums(fj) + (nvar/fixed_tau)/sample_size

  genet_var = tot_by_trait - error_var

  pve_j <- sj/genet_var

  return(list(h2 = genet_var, error_var = error_var, pve = pve_j))
}




