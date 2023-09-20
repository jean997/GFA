#'@export
gfa_estL_posthoc <- function(Y, fit, tol = 1e-5){
  flash_fit <- fit$fit$flash_fit
  n_new <- nrow(Y)
  s <- 1/sqrt(flash_fit_get_tau(flash_fit))

  nfct <- fit$fit$n_factors
  Lrand <- matrix(rnorm(n = n_new*nfct), nrow = n_new)
  g_ebnm <- list()
  fit_new <- flash_init(data = Y, S = s, var_type = NULL)
  for(i in seq(nfct)){
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
    flash_factors_fix(kset = seq(nfct), which_dim = "factors") %>%
    flash_backfit(tol = tol)
  if(length(fit$fixed_ix) > 0){
    L_hat <- fit_new$L_pm[,-fit$fixed_ix ]
    L_lfsr <- fit_new$L_lfsr[, -fit$fixed_ix]
  }else{
    L_hat <- fit_new$L_pm
    L_lfsr <- fit_new$L_lfsr
  }
  return(list(L_hat = L_hat, L_lfsr = L_lfsr))
}
