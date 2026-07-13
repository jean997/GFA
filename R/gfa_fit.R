
#'@title Genetic Factor Analysis
#'@param Z_hat A matrix of z-scores with rows for variants and columns for traits.
#'@param N Vector of sample sizes length equal to number of traits. If not provided,
#'N will default to the vector of 1's and the factor matrix will be returned on the "z-score scale".
#'@param N_case If all traits are continuous, omit this option. If some traits are binary, N_case should
#'be a vector with length equal to number of traits. Values should be NA for continuous traits or the number
#'of cases for binary traits.
#'@param pop_prev If all traits are continuous, omit this option. If some traits are binary, pop_prev should
#'be a vector with length equal to number of traits. Values should be NA for continuous traits or the population
#'prevalence for binary traits.
#'@param B_hat A matrix of GWAS effect estimates. B_hat is an alternative to Z_hat (only provide one of these).
#'If using B_hat, you must also provide S.
#'@param S If using B_hat, provide the corresponding matrix of standard errors.
#'@param R Estimated residual correlation matrix. This can be produced for example using R_ldsc or R_pt.
#'@param params List of parameters. For most users this can be left at the default values. See Details.
#'@param no_wrapup If TRUE, GFA will not perform wrap-up steps. Advanced option used for debugging.
#'@param F_init Initial estimate of F (optional).
#'@param fix_F,freeze_F Options for fixing F at initialization. See Details.
#'@return A list with elements L_hat and F_hat for estimated variant-factor and factor-trait effects, gfa_pve
#'which contains the proportion of heritability explained by each factor, and some other objects useful internally.
#'@details
#'You can fit GFA using two types of data. Either supply Z_hat and N or supply B_hat and S.
#'In either case, GFA will be fit with z-scores as inputs and factors will be scaled to
#'the standardized trait scale or the liability scale for binary traits. If you supply B_hat
#'and S, GFA will guess the scaling factors by computing the median of the ratio of each
#'trait's standard errors compared with the first trait's standard errors.
#'If you used B_hat/S and there are some binary traits, you only need to supply pop_prev.
#'If you use Z_hat/N and there are some binary traits, you must supply both N_case and pop_prev.
#'We recommend using the Z_hat/N option if sample sizes are available.
#'
#'The params list includes the following elements which can be modified.
#'Most users will not need to modify any of these, except possibly `max_iter`.
#'
#'kmax: Maximum number of factors. Defaults to twice the number of traits
#'
#'cond_num: Maximum allowable condition number for R. Defaults to 1000.
#'
#'max_iter: Maximum number of iterations. Defaults to 1000.
#'
#'extrapolate: Passed to flashier, defaults to TRUE
#'
#'ebnm_fn_F and ebnm_fn_L: Prior families for factors and loadings. Defaults to point-normal.
#'
#'init_fn: Flashier initialization function.
#'
#'fix_F and freeze_F provide different ways to fix initial estimates. If fix_F is TRUE, the factors
#'supplied to F_init will not be updated. However, GFA will be allowed to add additional factors to the final fit.
#'If freeze_F is TRUE, GFA will not add additional factors and the final value of F_hat will be equal to F_init
#'up to column scaling constants.
#'
#'The returned object will be a list with the following elements: F_hat (factor estimate), L_hat (loadings estimate), F_hat_single (the columns of
#'F_hat corresponding to single trait factors), F_hat_multi (the columns of F corresponding to multi-trait factors), fit (a flashier object),
#'scale (the scaling factor used), method (internal method type). If there are more than zero factors, the object will also include gfa_pve which
#'includes genet_var (the total variance explained by the set of variants used in estimation) and pve which is a traits by factors matrix.
#'The (i,j) element of pve gives the proportion of trait i hertiability explained by factor j.
#'
#'To compute credible intervals for pve, see gfa_intervals().
#'To compute GLS estimates of factor loadings see gfa_loadings_gls().
#'
#'@export
gfa_fit <- function(Z_hat = NULL,
                    N = NULL,
                    N_case = NULL,
                    pop_prev = NULL,
                    B_hat = NULL,
                    S = NULL,
                    R = NULL,
                    params = gfa_default_parameters(),
                    no_wrapup = FALSE,
                    F_init = NULL,
                    fix_F = FALSE,
                    freeze_F = FALSE){


  ## process parameters
  default_params <- gfa_default_parameters()
  method = "fixed_factors"
  #mode = "z-score"

  for(n in names(default_params)){
    if(is.null(params[[n]])) params[[n]] <- default_params[[n]]
  }
  for(n in names(params)){
    if(! n %in% names(default_params)) stop("Unknown parameter ", n, " provided.")
  }

  ## process inputs
  check_args_znbs(is.null(Z_hat), is.null(N), is.null(B_hat), is.null(S))
  if(is.null(Z_hat)){
    Z_hat <- B_hat/S
    scale <- get_scale_from_S(S)
    if(!is.null(pop_prev)){
      if(!length(pop_prev) == length(scale)){
        stop(paste0("pop_prev did not have expected length ", length(scale)))
      }
      ix <- which(!is.na(pop_prev))
      if(length(ix) > 0) scale[ix] <- scale[ix]/dnorm(qnorm(pop_prev[ix], lower.tail = TRUE))
    }
  }else{
    p <- ncol(Z_hat)
    if(is.null(N)){
      warning("Sample size not provided. Factor effects will be on the z-score scale which is sensitive to sample size.")
      N <- rep(1, ncol(Z_hat))
      scale <- sqrt(N)
    }else if(!is.null(N_case)){
      if(is.null(pop_prev)){
        stop("If supplying N_case, please also supply pop_prev")
      }
      if(!length(N) == p  & length(N_case) == p & length(pop_prev) == p){
        stop(paste0("N, N_case, and pop_prev do not all have expected length ", p))
      }
      scale <- sqrt(N)
      ix <- which(!is.na(N_case))
      if(length(ix) > 0) scale[ix] <- binary_const(N = N[ix], N_case = N_case[ix], pop_prev = pop_prev[ix])
    }else{
      if(!length(N) == p ){
        stop(paste0("N does not have expected length ", p))
      }
      scale <- sqrt(N)
    }
  }
  if((fix_F | freeze_F) & is.null(F_init)){
    stop("Cannot fix/freeze F, F_init not supplied.")
  }

  #dat <- gfa_set_data(Y = Z_hat, scale = scale, R = R, params = params, mode = mode)
  dat <- gfa_set_data(Y = Z_hat, scale = scale, R = R, params = params)


  dat$F_init <- F_init
  dat$fix_F <- fix_F
  dat$freeze_F <- freeze_F
  if(dat$freeze_F) dat$fix_F <- TRUE



  #method <- match.arg(method)

  if(is.null(dat$params$kmax)) dat$params$kmax <- 2*dat$p


  if(is.null(dat$R)){
    message("R is not supplied, fitting assuming full independence")
    fit <- fit_gfa_noR(dat)
    method <- "noR"
  }else if(method == "fixed_factors"){
    fit <- fit_gfa_ff(dat)
  }
  # }else if(method == "random_effect"){
  #   fit <- fit_gfa_re(dat)
  # }

  fit$method <- method
  ## wrap up
  if(is.null(fit$flash_fit$maxiter.reached) & !no_wrapup){

    fit <- fit %>% flash_nullcheck(remove = TRUE)
    fit <- gfa_duplicate_check(fit,
                               dim = 2,
                               check_thresh = params$duplicate_check_thresh)

    ret <- gfa_wrapup(fit,
                      method = method,
                      scale = dat$scale,
                      nullcheck = TRUE)
    #ret$mode <- mode
    ret$R <- dat$R
    ret$params <- dat$params
  }else{
    ret <- list(fit = fit,
                params = dat$params,
                scale = dat$scale,
                #mode = mode,
                R = dat$R)
  }
  return(ret)
}

fit_gfa_ff <- function(dat){
  eS <- eigen(dat$R)
  lambda_min <- eS$values[dat$p]
  vals <- eS$values - lambda_min
  if (dat$params$max_lr_percent < 1) {
    vv <- cumsum(vals)/sum(vals)
    nmax <- min(which(vv > dat$params$max_lr_percent))
  }else{
    nmax <- dat$p - 1
  }
  W <-eS$vectors[, seq(nmax), drop = FALSE] %*% diag(sqrt(vals[seq(nmax)]), ncol = nmax)

  # Initialize data object
  fit <-  flash_init(data = dat$Y, S = sqrt(lambda_min), var_type = dat$params$var_type)

  ## Initialize non-fixed factors.
  fit <- gfa_init_FL(fit, dat, noR = FALSE)
  k <- fit$n_factors

  ## Add fixed factors for R
  # randomly initialize A
  fixed_ebnm = flash_ebnm(prior_family = "normal",
                          g_init = ashr::normalmix(pi = c(1),mean = c(0),sd = c(1)),
                          fix_g = TRUE)
  A_rand <- matrix(rnorm(n = dat$n * nmax, mean = 0, sd = 1), nrow = dat$n, ncol = nmax)

  fit <- fit %>%
    flash_factors_init(., init = list(A_rand, W), ebnm_fn = fixed_ebnm) %>%
    flash_factors_fix(., kset = k + (1:nmax), which_dim = "factors") %>%
    flash_backfit(maxiter = dat$params$max_iter,
                  extrapolate=dat$params$extrapolate)
  fit$method <- "fixed_factors"
  return(fit)

}

fit_gfa_noR <- function(dat){
  #First initialize flash objects
  fit <-  flash_init(data = dat$Y, S = dat$S, var_type = dat$params$var_type)

  ## First initialize non-fixed factors.
  fit <- gfa_init_FL(fit, dat, noR = TRUE)
  k <- fit$n_factors

  # Backfit
  fit <- fit %>%
    flash_backfit(maxiter = dat$params$max_iter, extrapolate = dat$params$extrapolate)
  fit$method <- "noR"
  return(fit)
}


gfa_init_FL <- function(fit, dat, noR = FALSE){
  if(noR){
    Rinv <- diag(dat$p)
  }else{
    Rinv <- solve(dat$R)
  }
  if(!is.null(dat$F_init)){
    if(ncol(dat$F_init) <= ncol(dat$Y)){
      ftf_inv <- solve(t(dat$F_init) %*%  Rinv %*% dat$F_init)
    }else{
      ftf_inv <- corpcor::pseudoinverse(t(dat$F_init) %*% Rinv %*% dat$F_init)
    }

    ftx <- t(dat$F_init) %*% Rinv %*% t(dat$Y)
    L_init <- t(ftf_inv %*% ftx)

    fit <- fit %>%
      flash_factors_init(., init = list(L_init, dat$F_init),
                         ebnm_fn = list(dat$params$ebnm_fn_L, dat$params$ebnm_fn_F))
    if(dat$fix_F){
      fit <- fit %>%
        flash_factors_fix(., kset = 1:ncol(dat$F_init), which_dim = "factors")
    }
  }
  if(!dat$freeze_F){
    ## Add greedy factors if not frozen
    fit <-  fit %>%
      flash_greedy(Kmax = dat$params$kmax,
                   init_fn = dat$params$init_fn,
                   ebnm_fn = list(dat$params$ebnm_fn_L, dat$params$ebnm_fn_F) )
  }
  return(fit)
}

