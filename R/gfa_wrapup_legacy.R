#'@export
gfa_wrapup_legacy <- function(fit,  nullcheck = TRUE){
  if(nullcheck){
    fit <- fit %>% flash_nullcheck(remove = TRUE)
  }
  if(length(fit$flash_fit$fix.dim) == 0){
    fixed_ix <- rep(FALSE, ncol(fit$F_pm))
  }else{
    fixed_ix <- fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)})
  }

  if(any(fixed_ix)){
    fixed_ix <- which(fixed_ix)
    F_hat <- fit$F_pm[,-fixed_ix, drop=FALSE]
    L_hat <- fit$L_pm[, -fixed_ix, drop=FALSE]
  }else{
    F_hat <- fit$F_pm
    L_hat <- fit$L_pm
  }
#
#
#   Yhat <- L_hat %*% t(F_hat)
#   Lx <- colSums(L_hat^2)
#   Fx <- colSums(F_hat^2)
#   F_hat_scaled <- t(t(F_hat)/sqrt(Fx))
#   L_hat_scaled <- t(t(L_hat)/sqrt(Lx))
#   scale <- sqrt(Fx)*sqrt(Lx)
  ret <- list(fit=fit, L_hat = L_hat, F_hat = F_hat)
  return(ret)
}

#'@export
gfa_rebackfit_legacy <- function(fit, params){
  fit <- fit %>% flash_backfit(maxiter = params$max_iter,
                               extrapolate = params$extrapolate)
  if(is.null(fit$flash_fit$maxiter.reached)){
    ret <- gfa_wrapup_legacy(fit, nullcheck = TRUE)
  }else{
    ret <- list(fit = fit)
  }
  ret$params <- params
  return(ret)
}




