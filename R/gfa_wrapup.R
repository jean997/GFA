#'@export
gfa_wrapup <- function(fit, method, scale = NULL, num_single_fixed = 0, nullcheck = FALSE){
  F_hat_est <- fit$F_pm
  L_hat_est <- fit$L_pm

  if(!is.null(scale)){
    F_hat_est <- F_hat_est/scale
  }

  row_scale <- sqrt(colSums(F_hat_est^2))
  F_hat_est <- t(t(F_hat_est)/row_scale)
  L_hat_est <- t(t(L_hat_est)*row_scale)
  nfactor <- ncol(F_hat_est)

  fix.dim <- flashier:::get.fix.dim(fit$flash_fit)
  if(length(fix.dim) == 0){
    fixed_ix <- rep(FALSE, nfactor)
  }else{
    fixed_ix <- sapply(fix.dim, function(x){
      if(is.null(x)) return(FALSE)
      if(x == 2) return(TRUE)
      return(FALSE)})
  }

  est_ix <- which(!fixed_ix)
  if(nullcheck){
    fit <- fit %>% flash_nullcheck(tol = 0, remove = FALSE) # remove = FALSE to save indices
  }

  if(any(fixed_ix)){
    if(num_single_fixed > 0){
      single_ix <- (max(est_ix) + 1):(max(est_ix) + num_single_fixed)
      if(max(single_ix) < nfactor){
        error_ix <- (max(single_ix) + 1):length(fixed_ix)
      }else{
        error_ix <- NULL
      }
    }else{
      error_ix <- (max(est_ix) + 1):length(fixed_ix)
      single_ix <- NULL
    }
  }else{
    single_ix <- NULL
    error_ix <- NULL
  }

  if(any(fit$flash_fit$is.zero)){
    est_ix <- est_ix[!fit$flash_fit$is.zero[est_ix]]
    n <- length(est_ix)
    if(!is.null(single_ix)){
      single_ix <- single_ix[!fit$flash_fit$is.zero[single_ix]]
      n <- n + length(single_ix)
    }
    if(!is.null(error_ix)){
      error_ix <- error_ix[!fit$flash_fit$is.zero[error_ix]]
    }
    fit <- flash_factors_remove(fit, kset = which(fit$flash_fit$is.zero))
    error_ix <- (n+1):ncol(fit$F_pm) # update error_ix after removing null factors
  }

  F_hat <- F_hat_est[, est_ix, drop = FALSE]
  L_hat <- L_hat_est[, est_ix, drop = FALSE]

  F_hat_single <- NULL
  if(!is.null(single_ix)){
    F_hat_single <- F_hat_est[, single_ix, drop = FALSE]
  }

  ret <- list(fit=fit,
              method = method,
              L_hat = L_hat,
              F_hat = F_hat,
              F_hat_single = F_hat_single,
              num_single = length(single_ix),
              error_ix = error_ix,
              scale = scale)
  if(ncol(F_hat) > 0){
    ret$gfa_pve <- pve2(ret, error_ix)
  }
  return(ret)
}

#'@export
gfa_rebackfit <- function(gfa_fit, params){
  method <- gfa_fit$method
  scale <- gfa_fit$scale
  fit <- gfa_fit$fit %>% flash_backfit(maxiter = params$max_iter,
                               extrapolate = params$extrapolate)
  fit$method <- method
  if(is.null(fit$flash_fit$maxiter.reached)){
    fit <- fit %>% flash_nullcheck(remove = TRUE) #, tol = -Inf) # this will only remove 0 factors
    fit <- gfa_duplicate_check(fit,
                               dim = 2,
                               check_thresh = params$duplicate_check_thresh)
    ret <- gfa_wrapup(fit, method = method,
                      scale = dat$scale, nullcheck = FALSE)
    ret$params <- params
    #ret$mode <- fit$mode
    ret$R <- fit$R
  }else{
    ret <- list(fit = fit,
                params = fit$params, scale = fit$scale,
                #mode = fit$mode,
                R = fit$R)
  }
  return(ret)
}

### Scale notes:
## Standardized Effect model
## Bhatstd = L F^T + E
## Zhat = Bhatstd*diag(sqrt(N)) = L ( diag(sqrt(N)) F)^T + E
## So we have estimated F_hat_est = diag(sqrt(N)) F
## To retrieve F, F_hat = F_hat_est/sqrt(N)
## scale = sqrt(N)


