#'@title GFA Credible Intervals
#'@description Compute credible intervals for various estimated parameters.
#'@param gfa_fit An object output by gfa_fit.
#'@param nsamp Number of posterior samples.
#'@param level Level for credible intervals (e.g. 0.05 results in 95% credible intervals).
#'@param type A vector including some or all elements "pve", "F", or "L", see details.
#'@param variant_ix Optional list of indices if type includes "L".
#'@details
#'This function will compute credible intervals for the PVE matrix, F, or L depending on the
#'choice in "type". The returned object will include <type>_lower and <type>_upper for
#'each element in type.
#'
#'Please note that computing credible intervals for all elements of L can be very
#'time consuming, since L is typically very large. Your needs may be better served
#'by computing frequentist estimates and p-values for L using gfa_loadings_gls().
#'
#'@export
gfa_credints <- function(gfa_fit,
                         nsamp = 500,
                         level = 0.05,
                         type = c( "pve","F", "L"),
                         variant_ix = NULL){
  samps <- gfa_fit$fit$sampler(n = nsamp)
  plower <- level/2
  pupper <- 1-plower

  n <- nrow(samps[[1]]$L)
  p <- nrow(samps[[1]]$F)
  k <- ncol(samps[[1]]$L)

  F_orig <- gfa_fit$fit$F_pm
  col_scale <- gfa_fit$scale
  row_scale <- sqrt(colSums((F_orig/col_scale)^2))

  samps <- lapply(samps, function(s){
    s$L <- t(t(s$L)*row_scale)
    s$F <- t(t(s$F/col_scale)/row_scale)
    return(s)
  })

  if(gfa_fit$method == "fixed_factors"){
    fixed_ix <- which(gfa_fit$fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    est_ix <- which(!gfa_fit$fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
  }else{
    est_ix <- 1:k
  }

  ret <- list()
  if("F" %in% type){
    message("Computing intervals for F\n")

    F_ci <- purrr::map(samps, "F") %>%
            simplify2array() %>%
            apply(., 1:2, quantile, prob = c(plower, pupper))

    ret$F_lower <- F_ci[1,, est_ix]
    ret$F_upper <- F_ci[2,,est_ix]
  }
  if("L" %in% type){
    message("Computing intervals for L\n")
    warning("Computing credible intervals for L may take some time if L is large.\n")
    if(!is.null(variant_ix)){
      Ls <- purrr::map(samps, function(s){
        s$L[variant_ix,]
      })
      L_ci <- Ls %>%
        simplify2array() %>%
        apply(., 1:2, quantile, prob = c(plower, pupper))
    }else{
      L_ci <- purrr::map(samps, "L") %>%
              simplify2array() %>%
              apply(., 1:2, quantile, prob = c(plower, pupper))
    }
    ret$L_lower <- L_ci[1,, est_ix]
    ret$L_upper <- L_ci[2,,est_ix]
  }
  if("pve" %in% type){
    message("Computing intervals for pve\n")
    pve <- lapply(samps, function(s){
      pve_from_sample(gfa_fit, s)
    })
    pve_ci <- purrr::map(pve, "pve") %>%
              simplify2array() %>%
              apply(., 1:2, quantile, prob = c(plower, pupper))
    ret$pve_lower <- pve_ci[1,,]
    ret$pve_upper <- pve_ci[2,,]
  }
  return(ret)
}

pve_from_sample <- function(gfa_fit, samp){

  ef2 <- flashier:::get.EF2(gfa_fit$fit$flash_fit)
  k <- ncol(samp[[1]])
  p <- nrow(samp[[2]])
  n <- nrow(samp[[1]])

  if(gfa_fit$method == "fixed_factors"){
    fixed_ix <- which(gfa_fit$fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    est_ix <- which(!gfa_fit$fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    #cat("est_ix:", est_ix, "\n")
    # relative variance explained per trait/factor
    sj <- sapply(est_ix, function(kk){
      Vk <- flashier:::lowrank.expand(list((samp[[1]][,kk,drop = F])^2, (samp[[2]][,kk, drop = F])^2)) ## total effects on z-score scale
      ## standardized effect scale total variance
      colSums(Vk) # /(gfa_fit$scale^2) #sample should already be scaled
    }) |> matrix(nrow = p, byrow = FALSE)

    # variance from fixed factors (this is part of error)
    if(length(fixed_ix) == 0){
      fj_sums <- rep(0, p)
    }else{
      fj_sums <- sapply(fixed_ix, function(kk){
        Vk <- flashier:::lowrank.expand(list(ef2[[1]][,kk,drop = F], ef2[[2]][,kk, drop = F]))
        ## standardized effect scale total variance
        colSums(Vk)/(gfa_fit$scale^2)
      })|> matrix(nrow = p, byrow = FALSE) |> rowSums()
    }

    tau <- flashier:::get.tau(gfa_fit$fit$flash_fit)
    if(inherits(tau, "matrix")){
      stopifnot(nrow(tau) == n & ncol(tau) == p) # should never fail
      S2 <- gfa_fit$fit$flash_fit$given.S2
      est_tau <- 1/((1/tau) - S2)
    }else{
      stopifnot(length(tau) == p) # this should never fail
      # error component of tau
      fixed_tau <- flash_fit_get_fixed_tau(gfa_fit$fit$flash_fit)
      est_tau <- 1/((1/tau) - (1/fixed_tau))
    }
    error_var <- fj_sums + (n/fixed_tau)/(gfa_fit$scale^2)
    total_var <- rowSums(sj) + fj_sums + (n/tau)/(gfa_fit$scale^2)

    genet_var = rowSums(sj) + (n/est_tau)/(gfa_fit$scale^2)
    pve_j <- sj/genet_var
  }else if(gfa_fit$method == "noR"){
    est_ix <- 1:k
    # relative variance explained per trait/factor
    sj <- sapply(est_ix, function(kk){
      Vk <- flashier:::lowrank.expand(list(samp[[1]][,kk,drop = F]^2, samp[[2]][,kk, drop = F]^2)) ## total effects on z-score scale
      ## standardized effect scale total variance
      colSums(Vk) #/(gfa_fit$scale^2)
    }) |> matrix(nrow = p, byrow = FALSE)

    tau <- flashier:::get.tau(gfa_fit$fit$flash_fit)
    if(inherits(tau, "matrix")){
      stopifnot(nrow(tau) == n & ncol(tau) == p) # should never fail
      S2 <- gfa_fit$fit$flash_fit$given.S2
      error_var <-  colSums(S2)/(gfa_fit$scale^2)
      tau_var <- colSums(1/tau)
      gvar_tau <- tau_var - error_var
    }else{
      stopifnot(length(tau) == p) # this should never fail
      # error component of tau
      fixed_tau <- flash_fit_get_fixed_tau(gfa_fit$fit$flash_fit)
      est_tau <- 1/((1/tau) - (1/fixed_tau))
      error_var <-  (n/fixed_tau)/(gfa_fit$scale^2)
      gvar_tau <- (n/est_tau)/(gfa_fit$scale^2)
    }
    #total_var <- rowSums(sj) + (n/tau)/(gfa_fit$scale^2)

    genet_var = rowSums(sj) + gvar_tau
    pve_j <- sj/genet_var
  }
  return(list(genet_var = genet_var, pve = pve_j))
}
