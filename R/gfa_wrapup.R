#'@export
gfa_wrapup <- function(fit, method, scale = NULL,
                       n_single = 0,
                       nullcheck = FALSE){

  fix_ix <- my_flash_get_fixed_idx(fit)
  if(!is.null(fix_ix$loadings)){
    stop("Something is wrong. Loadings are fixed.")
  }
  ntotal <- fit$n_factors
  n_est <- ntotal - length(fix_ix$factors)
  n_error <- ntotal - n_est - n_single

  if(nullcheck){
    fit <- gfa_nullcheck(fit, num_single_fixed = n_single, num_error_fixed =n_error)
    n_single <- fit$n_single
    n_est <- fit$n_est
    n_error <- fit$n_error
    ntotal <- n_est + n_single + n_error
    fit <- fit$fit
  }


  F_hat_est <- fit$F_pm
  L_hat_est <- fit$L_pm

  if(!is.null(scale)){
    F_hat_est <- F_hat_est/scale
  }

  row_scale <- sqrt(colSums(F_hat_est^2))
  F_hat_est <- t(t(F_hat_est)/row_scale)
  L_hat_est <- t(t(L_hat_est)*row_scale)


  if(n_est > 0){
    est_ix <- 1:n_est
    F_hat <- F_hat_est[, est_ix, drop = FALSE]
    L_hat <- L_hat_est[, est_ix, drop = FALSE]
  }else{
    F_hat <- NULL
    L_hat <- NULL
  }
  if(n_single > 0){
    single_ix <- n_est + (1:n_single)
    F_hat_single <- F_hat_est[, single_ix, drop = FALSE]
  }else{
    single_ix <- c()
    F_hat_single <- NULL
  }
  if(n_error > 0){
    error_ix <- n_est + n_single + (1:n_error)
  }else{
    error_ix <- c()
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
gfa_rebackfit <- function(gfa_fit, params, single_check = TRUE, wrapup = TRUE){
  method <- gfa_fit$fit$method
  scale <- gfa_fit$scale
  fit <- gfa_fit$fit %>% flash_backfit(maxiter = params$max_iter,
                               extrapolate = params$extrapolate)
  fit$method <- method

  if(is.null(fit$flash_fit$maxiter.reached) & wrapup){
    fit <- flash_nullcheck(fit, tol = 0, remove = TRUE)
    fit <- gfa_duplicate_check(fit,
                               dim = 2,
                               check_thresh = params$duplicate_check_thresh)
    if(single_check){
      fit <- gfa_singletrait_check(fit, check_thresh = params$singletrait_check_thresh, params = params)

      ret <- gfa_wrapup(fit,
                        method = method,
                        scale = scale,
                        num_single_fixed = fit$num_single_fixed,
                        nullcheck = TRUE)
    }else{
      ret <- gfa_wrapup(fit,
                        method = method,
                        scale = scale,
                        num_single_fixed = 0,
                        nullcheck = TRUE)
    }
    ret$R <- gfa_fit$R
    ret$params <- params
  }else{
    ret <- list(fit = fit,
                params = dat$params,
                scale = dat$scale,
                R = dat$R)
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


