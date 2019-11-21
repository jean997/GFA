
mask_flashier <- function(n, mask_proportion, mats, var_type  =  c("constant", "by_row", "by_col", "kronecker",
                                               "zero", "noisy_constant", "noisy_byrow",
                                               "noisy_bycol") ,
                          init_type = c("flashier", "soft_impute", "from_data")){
  var_type <- match.arg(var_type)
  init_type <- match.arg(init_type)
  mats$beta_hat[is.na(mats$se_hat)] <- NA
  args <- list(data = mats$beta_hat, output.lvl = 1, tol = 0.01)


  fit_func <- function(mats){
    f <- run_flashier(mats, var_type, init_type)
    fact <- f$loadings.pm[[2]]
    load <- f$loadings.pm[[1]]
    D <- f$loadings.scale
    EY <- load %*% diag(D) %*% t(fact)
    return(EY)
  }

  mask_res <- mask_and_fit(mats, fit_func, mask_proportion=mask_proportion, n_times=n)
}
