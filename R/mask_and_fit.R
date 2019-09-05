#'@title Mask data and Fit
#'@param mats List containing beta_hat and se_hat matrices
#'@param fit_function Function that takes mats as input and returns E[beta_hat]
#'@param mask_proportion Proportion of non-missing values to mask
#'@param n_times Number of replications
#'@param seed seed
#'@export
mask_and_fit <- function(mats, fit_function,
                         mask_proportion= 0.05, n_times = 5, seed = 100){

  set.seed(seed)
  non_missing_index <- which(!is.na(mats$beta_hat))
  n_mask <- round(mask_proportion*length(non_missing_index))
  mask_ix <- replicate(n=n_times, {sample(non_missing_index, size=n_mask)})
  errs <- apply(mask_ix, 2, function(ix){
                      my_mats <- mats
                      my_mats$beta_hat[ix] <- NA
                      my_mats$se_hat[ix] <- NA
                      est <- fit_function(my_mats)
                      sum((mats$beta_hat[ix]-est[ix])^2)
                    })
  return(list(errs = errs, mask_ix = mask_ix))

}

#'@export
fit_flash_zero_svd <- function(mats){
  mats$beta_hat[is.na(mats$se_hat)] <- NA
  mats$se_hat[is.na(mats$beta_hat)] <- Inf
  kmax <- 100
  data <- with(mats, flash_set_data(Y = beta_hat, S = se_hat))
  f  <- flash_add_factors_from_data(data,K=kmax, var_type="zero")
  f <- flash_backfit(data,f, var_type="zero")
  D <- f$ldf$d
  fact <- f$ldf$f
  load <- f$ldf$l
  EY <- load %*% diag(D) %*% t(fact)
  return(EY)
}

#'@export
fit_flash_zero_greedy <- function(mats){
  mats$beta_hat[is.na(mats$se_hat)] <- NA
  mats$se_hat[is.na(mats$beta_hat)] <- Inf
  kmax <- 100
  data <- with(mats, flash_set_data(Y = beta_hat, S = se_hat))
  f  <- flash_add_greedy(data,K=kmax, var_type="zero")
  f <- flash_backfit(data,f, var_type="zero")
  D <- f$ldf$d
  fact <- f$ldf$f
  load <- f$ldf$l
  EY <- load %*% diag(D) %*% t(fact)
  return(EY)
}

#'@export
fit_flashier_null  <- function(mats){
  mats$beta_hat[is.na(mats$se_hat)] <- NA
  mats$se_hat[is.na(mats$beta_hat)] <- Inf
  kmax <- 100
  f <- with(mats, flashier(data=beta_hat, S = se_hat, var.type=NULL,
                             greedy.Kmax = kmax, init.tol=1e-4,
                             fit = "full"))

  fact <- f$loadings.pm[[2]]
  load <- f$loadings.pm[[1]]
  D <- f$loadings.scale
  EY <- load %*% diag(D) %*% t(fact)
  return(EY)
}

#'@export
fit_flashier_zero  <- function(mats){
  mats$beta_hat[is.na(mats$se_hat)] <- NA
  mats$se_hat[is.na(mats$beta_hat)] <- Inf
  kmax <- 100
  f <- with(mats, flashier(data=beta_hat, S = se_hat, var.type=0,
                           greedy.Kmax = kmax, init.tol=1e-4,
                           fit = "full"))

  fact <- f$loadings.pm[[2]]
  load <- f$loadings.pm[[1]]
  D <- f$loadings.scale
  EY <- load %*% diag(D) %*% t(fact)
  return(EY)
}

