#'@export
gfa_wrapup <- function(fit, method, scale = NULL, nullcheck = FALSE){
  if(nullcheck){
    fit <- fit %>% flash_nullcheck(remove = TRUE)
  }

  F_hat <- fit$F_pm
  L_hat <- fit$L_pm

  if(method == "fixed_factors"){
    fixed_ix <- which(fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    if(length(fixed_ix) > 0){
      F_hat <- fit$F_pm[,-fixed_ix, drop=FALSE]
      L_hat <- fit$L_pm[, -fixed_ix, drop=FALSE]
    }
  }
  #F_hat_est <- F_hat
  #L_hat_est <- L_hat

  if(!is.null(scale)){
    F_hat <- F_hat/scale
  }
  row_scale <- sqrt(colSums(F_hat^2))
  F_hat <- t(t(F_hat)/row_scale)
  L_hat <- t(t(L_hat)*row_scale)
  #L_hat_est <- t(t(L_hat_est)*row_scale)
  #F_hat_est <- t(t(F_hat_est)/row_scale)

  F_hat_multi <- F_hat
  F_hat_single <- NULL
  if(ncol(F_hat) > 0){
    hat_single <- find_single_trait(F_hat)
    if(length(hat_single) > 0){
      F_hat_single <- F_hat[,hat_single,drop = F]
      F_hat_multi <- F_hat[, -hat_single, drop = F]
    }
  }
  ret <- list(fit=fit,
              method = method,
              L_hat = L_hat,
              F_hat = F_hat,
              F_hat_single = F_hat_single,
              F_hat_multi = F_hat_multi,
              scale = scale)
  if(ncol(F_hat) > 0){
    ret$gfa_pve <- pve2(ret)
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


