## assumes order of estimated, single fixed, error fixed
gfa_nullcheck <- function(fit, num_single_fixed, num_error_fixed){
  stopifnot(inherits(fit, "flash"))
  n_est <- fit$n_factors - num_single_fixed - num_error_fixed
  if(length(n_est) > 0){
    est_ix <- 1:n_est
    n <- n_est
  }else{
    est_ix <- c()
    n <- 0
  }
  if(length(num_single_fixed) > 0){
    single_ix <- n + (1:num_single_fixed)
    n <- n + num_single_fixed
  }else{
    single_ix <- c()
  }
  if(length(num_error_fixed) > 0){
    error_ix <- n + (1:num_error_fixed)
  }else{
    error_ix <- c()
  }
  fit <- flash_nullcheck(fit, remove = FALSE)
  zero_ix <- which(fit$flash_fit$is.zero)
  if(length(zero_ix) > 0){
    n_est <- n_est - sum(est_ix %in% zero_ix)
    n_single <- num_single_fixed - sum(single_ix %in% zero_ix)
    n_error <- num_error_fixed - sum(error_ix %in% zero_ix)
    fit <- flash_factors_remove(fit, kset = zero_ix)
  }else{
    n_single <- num_single_fixed
    n_error <- num_error_fixed
  }

  return(list(fit = fit, n_est = n_est, n_single = n_single, n_error = n_error))
}
