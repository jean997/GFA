#'@export
gfa_wrapup <- function(fit, method, scale = NULL, nullcheck = TRUE){
  if(nullcheck){
    fit <- fit %>% flash_nullcheck(remove = TRUE)
  }


  if(method == "noR" | method == "random_effect"){
    F_hat <- fit$F_pm
    L_hat <- fit$L_pm
  }else if(method == "fixed_factors"){
    fixed_ix <- which(fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    F_hat <- fit$F_pm[,-fixed_ix, drop=FALSE]
    L_hat <- fit$L_pm[, -fixed_ix, drop=FALSE]
  }

  F_hat_est <- F_hat
  if(!is.null(scale)){
    F_hat <- F_hat*scale
  }
  if(ncol(F_hat) > 0){
    hat_single <- find_single_trait(F_hat)
    if(length(hat_single) > 0){
      F_hat_single <- F_hat[,hat_single]
      F_hat_multi <- F_hat[, -hat_single]
    }else{
      F_hat_single <- NULL
      F_hat_multi <- F_hat
    }
  }
  ret <- list(fit=fit, method = method,
              L_hat = L_hat, F_hat = F_hat,
              F_hat_single = F_hat_single,
              F_hat_multi = F_hat_multi,
              scale = scale, F_hat_est = F_hat_est)
  return(ret)
}

#'@export
gfa_rebackfit <- function(fit, params){
  fit <- fit %>% flash_backfit(maxiter = params$max_iter,
                               extrapolate = params$extrapolate)
  if(is.null(fit$flash_fit$maxiter.reached)){
    ret <- gfa_wrapup(fit, nullcheck = TRUE)
  }else{
    ret <- list(fit = fit)
  }
  ret$params <- params
  return(ret)
}




