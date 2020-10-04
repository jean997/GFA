#'@title Fit plain vanilla FLASH to standardized effects
#'@param B_hat Matrix of (non-standardized) effect estimates SNPs by traits
#'@param S_hat Matrix of standard errors SNPs by traits
#'@param N Vector of sample sizes length equal to number of traits
#'@param kmax Maximum number of factors
#'@param init_fn init.fn argument for flashier
#'@param prior_family prior.family argument for flashier
#'@return A list with elements fit, B_hat, L_hat, F_hat
#'@export
fit_plain <- function(B_hat, S_hat, N, kmax=100,
                      init_fn =init.fn.softImpute,
                      prior_family = prior.point.normal(),
                      adjust=TRUE){

  n_var <- nrow(B_hat)
  n_trait <- ncol(B_hat)
  stopifnot(nrow(S_hat) == n_var & ncol(S_hat) == n_trait)

  if(!missing(N)) stopifnot(length(N) == n_trait)
  if(adjust & missing(N)) stop("To adjust please supply N.")



  if(adjust){
    B_tilde = t( (1/sqrt(N)) *t(B_hat/S_hat))
    S_tilde = t( (1/sqrt(N)) * t(matrix(1, nrow=n_var, ncol=n_trait)))
  }else{
    B_tilde = B_hat
    S_tilde = S_hat
  }


  fit <- flash.init(data=B_tilde, S = S_tilde,  var.type=2) %>%
    flash.add.greedy(Kmax = kmax, init.fn = init_fn, prior.family = prior_family) %>%
    flash.backfit() %>%
    flash.nullcheck()

  if(fit$n.factors > 0){
    F_hat <- fit$loadings.pm[[2]]
    L_hat <- fit$loadings.pm[[1]]
    B_hat <- fitted(fit)
  }else{
    B_hat <- F_hat  <- L_hat <- NULL;
  }
  ret <- list(fit=fit, B_hat = B_hat, L_hat = L_hat, F_hat = F_hat)
  return(ret)
}
