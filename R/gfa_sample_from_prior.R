gfa_one_sample_from_prior <- function(gfa_obj){
  n <- nrow(gfa_obj$fit$L_pm)
  p <- nrow(gfa_obj$fit$F_pm)
  k <- ncol(gfa_obj$fit$L_pm)
  Fmat <- matrix(nrow = p, ncol = k)
  Lmat <- matrix(nrow = n, ncol = k)

  ## step 1 fill in fixed columns
  # note currently fix_ix$loadings is always null for gfa objects
  fix_ix <- my_flash_get_fixed_idx(gfa_obj$fit)
  if(length(fix_ix$factors) > 0){
    Fmat[,fix_ix$factors] <- gfa_obj$fit$F_pm[, fix_ix$factors]
  }

  ## step 2 sample nonfixed factors and loadings
  for(i in 1:k){
    # sample loadings
    l_ghat <- gfa_obj$fit$L_ghat[[i]]
    myl <- GWASBrewer::rnormalmix(n = n, sd = l_ghat$sd , pi = l_ghat$pi, mu = l_ghat$mean)
    Lmat[,i] <- myl
    # sample factor
    if(! i %in% fix_ix$factors){
      f_ghat <- gfa_obj$fit$F_ghat[[i]]
      myf <- GWASBrewer::rnormalmix(n = p, sd = f_ghat$sd , pi = f_ghat$pi, mu = f_ghat$mean)
      Fmat[,i] <- myf
    }
  }

  ## compute separately LF^T + Theta and E
  if(length(gfa_obj$error_ix) > 0){
    eix <- gfa_obj$error_ix
    lft <- Lmat[, -eix] %*% t(Fmat[, -eix])
  }else{
    lft <- Lmat %*% t(Fmat)
  }

  s_error <- 1/sqrt(flash_fit_get_fixed_tau(gfa_obj$fit$flash_fit)) ## s for last e.v.
  s <- 1/sqrt(flash_fit_get_tau(gfa_obj$fit$flash_fit)) # total s for last e.v. plus theta
  s_theta <- sqrt(s^2 - s_error^2)
  stopifnot(length(s) == p)

  # sample theta
  theta <- Reduce(cbind, lapply(s_theta, function(mysd){ rnorm(n = n, sd = mysd, mean = 0)}))

  # sample etilde
  etilde <- Reduce(cbind, lapply(s_error, function(mysd){ rnorm(n = n, sd = mysd, mean = 0)}))

  Y <- Lmat %*% t(Fmat) + theta + etilde
  colnames(Y) <- NULL
  Ymean <- lft + theta
  colnames(Ymean) <- NULL
  return(list(Y = Y, Ymean = Ymean))
}
