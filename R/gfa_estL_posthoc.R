#'@export
gfa_estL_theta_posthoc <- function(Y, fit, tol = 1e-5){
  flash_fit <- fit$fit$flash_fit
  n_new <- nrow(Y)
  s_error <- 1/sqrt(flash_fit_get_fixed_tau(flash_fit)) ## s for last e.v.

  s <- 1/sqrt(flash_fit_get_tau(flash_fit)) # total s for last e.v. plus theta
  s_theta <- sqrt(s^2 - s_error^2)

  nfactor <- fit$fit$n_factors
  ntrait <- nrow(fit$fit$F_pm)

  ## first fit without theta separated
  fit_new <- flash_init(data = Y, S = s, var_type = NULL)

  EF <- flash_fit_get_pm(flash_fit, 2)
  gL <- flash_fit_get_g(flash_fit, 1) # g for columns doesn't matter because they will always be fixed
  for(i in seq(nfactor)){
    if(class(gL[[i]]) == "normalmix"){
      if(length(gL[[i]]$mean) == 1){
        g_ebnm = flash_ebnm(prior_family = "normal",
                                 g_init = gL[[i]],
                                 fix_g = TRUE)
      }else{
        g_ebnm = flash_ebnm(prior_family = "normal_scale_mixture",
                                 g_init = gL[[i]],
                                 fix_g = TRUE)
      }
      #Lrandi <- GWASBrewer::rnormalmix(n = n_new, pi = gL[[i]]$pi, sd = gL[[i]]$sd, mu = gL[[i]]$mean)
      #Lrandi <- matrix(Lrandi, ncol = 1)
      Lrandi <- matrix(0, ncol = 1, nrow = n_new)
    }else if(class(gL) == "unimix"){
      stop("This function is only implemented for normal and normal mixture priors.")
      #g_pfam = prior.unimodal.symmetric(scale = 1, g_init=gL, fix_g = TRUE)
    }
    EFi <- EF[,i, drop = FALSE]
    fit_new <- fit_new %>%
      flash_factors_init(init = list(Lrandi, EFi),  ebnm_fn = g_ebnm)
  }
  fit_new <- fit_new %>%
    flash_factors_fix(kset = seq(nfactor), which_dim = "factors") %>%
    flash_backfit(tol = tol)

  # Now fit again with theta separate
  fit_new2 <- flash_init(data = Y, S = s_error, var_type = NULL)
  EF1 <- flash_fit_get_pm(fit_new$flash_fit, 1)
  EF2 <- flash_fit_get_pm(fit_new$flash_fit, 2)
  gL <- flash_fit_get_g(fit_new$flash_fit, 1)
  for(i in 1:nfactor){
    if(class(gL[[i]]) == "normalmix"){
      if(length(gL[[i]]$mean) == 1){
        g_ebnm = flash_ebnm(prior_family = "normal",
                            g_init = gL[[i]],
                            fix_g = TRUE)
      }else{
        g_ebnm = flash_ebnm(prior_family = "normal_scale_mixture",
                            g_init = gL[[i]],
                            fix_g = TRUE)
      }
    }
    fit_new2 <- fit_new2 %>%
      flash_factors_init(init = list(EF1[,i,drop = FALSE], EF2[,i,drop = FALSE]),  ebnm_fn = g_ebnm)
  }
  # Add factors for theta
  gL1 <- ashr::normalmix(pi = 1, mean = 0, sd = 1)
  for(i in seq(ntrait)){
     g_ebnm = flash_ebnm(prior_family = "normal",
                              g_init = gL1,
                              fix_g = TRUE)
     EFi <- matrix(0, nrow = ntrait, ncol = 1)
     EFi[i] <- s_theta[i]
     #Lrandi <- matrix(rnorm(n= n_new), ncol = 1)
     fit_new2 <- fit_new2 %>%
       flash_factors_init(init = list(matrix(0,nrow = n_new, ncol = 1), EFi),  ebnm_fn = g_ebnm)
  }
  fit_new2 <- fit_new2 %>%
     flash_factors_fix(kset = seq(nfactor + ntrait), which_dim = "factors") %>%
     flash_backfit(tol = tol) #, kset = seq(ntrait) + nfactor) perhaps a bug in kset argument

  ## Scale
  F_hat_est <- fit_new2$F_pm
  L_hat_est <- fit_new2$L_pm

  #if(!is.null(fit$scale)){
  #  F_hat_est <- F_hat_est/fit$scale
  #}

  #row_scale <- sqrt(colSums(F_hat_est^2))
  #F_hat_est <- t(t(F_hat_est)/row_scale)
  #L_hat_est <- t(t(L_hat_est)*row_scale)
  ####

  n_est <- nfactor - fit$num_single - length(fit$error_ix)
  n_single <- fit$num_single
  n_error <- length(fit$error_ix)
  n_theta <- ntrait

  est_ix <- seq(n_est)
  L_multi <- L_hat_est[, est_ix, drop = FALSE]
  F_multi <- F_hat_est[,est_ix, drop = FALSE]
  H_multi <- L_multi %*% t(F_multi)
  if(n_single > 0){
    single_ix <- seq(n_single) + n_est
    L_single <- L_hat_est[, single_ix, drop = FALSE]
    F_single <- F_hat_est[, single_ix, drop = FALSE]
    H_single <- L_single %*% t(F_single)
    H <- H_multi + H_single
  }else{
    H_single <- matrix(0, nrow = nrow(Y), ncol = ntrait)
    H <- H_multi
  }
  if(n_error > 0){
    error_ix <- seq(n_error) + n_est + n_single
  }
  theta_ix <- seq(n_theta) + n_est + n_single + n_error
  L_theta <- L_hat_est[, theta_ix, drop = FALSE]
  F_theta <- F_hat_est[, theta_ix, drop = FALSE]
  theta <- L_theta %*% t(F_theta)

  Y_pm <- H + theta
  cat("Please note that posterior means are returned on the z-score scale. Use scale object to rescale to standardized effects.")
  return(list(theta_pm = theta, H_pm = H_multi, H_single_pm = H_single, Y_pm = Y_pm, scale = fit$scale))
}


#'@export
gfa_estL_posthoc <- function(Y, fit, tol = 1e-5){
  flash_fit <- fit$fit$flash_fit
  n_new <- nrow(Y)
  s <- 1/sqrt(flash_fit_get_tau(flash_fit))

  nfactor <- fit$fit$n_factors
  Lrand <- matrix(rnorm(n = n_new*nfactor), nrow = n_new)
  g_ebnm <- list()
  fit_new <- flash_init(data = Y, S = s, var_type = NULL)
  for(i in seq(nfactor)){
    # g for columns doesn't matter because they will always be fixed
    gL <- flash_fit_get_g(flash_fit, 1)[[i]]
    if(class(gL) == "normalmix"){
      if(length(gL$mean) == 1){
        g_ebnm[[i]] = flash_ebnm(prior_family = "normal",
                                 g_init = gL,
                                 fix_g = TRUE)
      }else{
        g_ebnm[[i]] = flash_ebnm(prior_family = "normal_scale_mixture",
                                 g_init = gL,
                                 fix_g = TRUE)
      }
    }else if(class(gL) == "unimix"){
      stop("Not implemented yet.")
      g_pfam = prior.unimodal.symmetric(scale = 1, g_init=gL, fix_g = TRUE)
    }
    EFk <- flash_fit_get_pm(flash_fit, 2)[,i, drop = FALSE]
    fit_new <- fit_new %>%
      flash_factors_init(init = list(Lrand[,i, drop = FALSE], EFk),  ebnm_fn = g_ebnm[[i]])
  }
  fit_new <- fit_new %>%
    flash_factors_fix(kset = seq(nfactor), which_dim = "factors") %>%
    flash_backfit(tol = tol)

  fit_new <- gfa_wrapup(fit_new, method = fit$method, scale = fit$scale)
  return(fit_new)
}
