#'@export
gfa_estL_theta_posthoc <- function(Y, fit, tol = 1e-5){
  flash_fit <- fit$fit$flash_fit
  n_new <- nrow(Y)
  s_error <- 1/sqrt(flash_fit_get_fixed_tau(flash_fit)) ## s for last e.v.


  s <- 1/sqrt(flash_fit_get_tau(flash_fit)) # total s for last e.v. plus theta
  s_theta <- sqrt(s^2 - s_error^2)

  V_F_hat <- fit$fit$F_psd^2

  nfactor <- fit$fit$n_factors
  ntrait <- nrow(fit$fit$F_pm)

  ## refit with priors and F fixed to get loadings for Y.
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

  ## Get posterior of theta using ebnm

  R <- Y - with(fit_new, L_pm %*% t(F_pm))
  #R <- b_hat - myL %*% t(myF)

  V_lft <- with(fit_new, (L_psd^2) %*% t(V_F_hat) + (L_psd^2) %*% t(F_pm^2) + (L_pm^2) %*% t(V_F_hat))
  tau_R <- sqrt(s_error^2 + V_lft)

  ebnm_fit <- lapply(seq(ntrait), function(j){
    g <- ashr::normalmix(pi = 1, mean = 0, sd = s_theta[j])
    ebnm(x = R[,j], s = tau_R[,j], prior_family = "normal",
         g_init = g, fix_g = TRUE,
         optmethod = "nlm")
  })
  theta_pm <- lapply(ebnm_fit, function(f){f$posterior$mean}) %>% reduce(., cbind)
  V_theta <- lapply(ebnm_fit, function(f){f$posterior$sd^2}) %>% reduce(., cbind)


  n_est <- nfactor - fit$num_single - length(fit$error_ix)
  n_single <- fit$num_single
  n_error <- length(fit$error_ix)

  ## dont scale
  F_hat_est <- fit_new$F_pm
  L_hat_est <- fit_new$L_pm

  V_L_hat <- fit_new$L_psd^2

  est_ix <- seq(n_est)
  L_multi <- L_hat_est[, est_ix, drop = FALSE]
  F_multi <- F_hat_est[,est_ix, drop = FALSE]
  H_multi <- L_multi %*% t(F_multi)

  V_L_multi <- V_L_hat[, est_ix, drop = FALSE]
  V_F_multi <- V_F_hat[, est_ix, drop = FALSE]
  V_H_multi <- V_L_multi %*% t(V_F_multi) + V_L_multi %*% t(F_multi^2) + (L_multi^2) %*% t(V_F_multi)

  if(n_single > 0){
    single_ix <- seq(n_single) + n_est
    L_single <- L_hat_est[, single_ix, drop = FALSE]
    F_single <- F_hat_est[, single_ix, drop = FALSE]
    H_single <- L_single %*% t(F_single)

    V_L_single <- V_L_hat[, single_ix, drop = FALSE]
    V_H_single <- V_L_single %*% t(F_single^2)

    H <- H_multi + H_single
    V_H <- V_H_multi + V_H_single
  }else{
    H_single <- V_H_single <- matrix(0, nrow = nrow(Y), ncol = ntrait)
    H <- H_multi
    V_H <- V_H_multi
  }

  Y_pm <- H + theta_pm
  V_Y <- V_H + V_theta
  cat("Please note that posterior means are returned on the z-score scale. Rescale by multiplying by standard errors.\n")
  return(list(H_pm = H_multi,
              H_single_pm = H_single,
              theta_pm = theta_pm,
              Y_pm = Y_pm,
              V_H_multi = V_H_multi,
              V_H_single = V_H_single,
              V_theta = V_theta,
              V_Y = V_Y))
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
