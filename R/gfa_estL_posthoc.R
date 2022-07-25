#'@export
gfa_estL_posthoc <- function(Y, fit, tol = 1e-5){
  flash_fit <- fit$fit$flash.fit
  n_new <- nrow(Y)
  s <- 1/sqrt(ff.tau(flash_fit))

  nfct <- fit$fit$n.factors
  Lrand <- matrix(rnorm(n = n_new*nfct), nrow = n_new)
  g_ebnm <- list()
  fit_new <- flash.init(data = Y, S = s, var.type = NULL)
  for(i in seq(nfct)){
    # g for columns doesn't matter because they will always be fixed
    gL <- ff.g(flash_fit, 1)[[i]]
    if(class(gL) == "normalmix"){
      if(length(gL$mean) == 1){
        g_ebnm[[i]] = as.ebnm.fn(prior_family = "normal",
                                 g_init = gL,
                                 fix_g = TRUE)
      }else{
        g_ebnm[[i]] = as.ebnm.fn(prior_family = "normal_scale_mixture",
                                 g_init = gL,
                                 fix_g = TRUE)
      }
    }else if(class(gL) == "unimix"){
      stop("Not implemented yet.")
      g_pfam = prior.unimodal.symmetric(scale = 1, g_init=gL, fix_g = TRUE)
    }
    EFk <- ff.pm(flash_fit, 2)[,i, drop = FALSE]
    fit_new <- fit_new %>%
      flash.init.factors(init = list(Lrand[,i, drop = FALSE], EFk),  ebnm.fn = g_ebnm[[i]])
  }
  fit_new <- fit_new %>%
    flash.fix.factors(kset = seq(nfct), mode = 2) %>%
    flash.backfit(tol = tol)
  if(length(fit$fixed_ix) > 0){
    L_hat <- fit_new$L.pm[,-fit$fixed_ix ]
    L_lfsr <- fit_new$L.lfsr[, -fit$fixed_ix]
  }else{
    L_hat <- fit_new$L.pm
    L_lfsr <- fit_new$L.lfsr
  }
  return(list(L_hat = L_hat, L_lfsr = L_lfsr))
}
