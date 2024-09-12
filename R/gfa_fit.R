
#'@title Genetic Factor Analysis
#'@param Z_hat A matrix of z-scores with rows for variants and columns for traits.
#'@param N Vector of sample sizes length equal to number of traits. If not provided,
#'N will default to the vector of 1's and the factor matrix will be returned on the "z-score scale".
#'@param B_hat A matrix of GWAS effect estimates. B_hat is an alternative to Z_hat (only provide one of these).
#'If using B_hat, you must also provide S.
#'@param S If using B_hat, provide the corresponding matrix of standard errors.
#'@param R Estimated residual correlation matrix. This can be produced for example using R_ldsc or R_pt.
#'@param params List of parameters. For most users this can be left at the default values. See details.
#'@param method Either "fixed_factors" or "random_effect". "random_effect" is experimental. See details.
#'@param mode Either "z-score" or "b-std". See details.
#'@return A list with elements L_hat and F_hat for estimated variant-factor and factor-trait effects, gfa_pve
#'which contains the proportion of heritability explained by each factor, and some other objects useful internally.
#'@details
#'
#'The method argument can be either fixed_factors or random_effect. These are different methods
#'for fitting the GFA model. The random_effect method requires installing a fork of the flashier
#'package using `install_github("jean997/flashier", ref = "random-effect")`. Most users can
#'stick with the fixed_factors method.
#'
#'The mode option tells GFA to either fit using z-scores as the outcome and
#'then convert the factor matrix to the standardized scale (mode = "z-score")
#'or to fit directly using standardized
#'effect estimates as the outcome (mode = "b-std"). These give
#'very similar results and the z-score mode is faster and so recommended.
#'
#'
#'The params list includes the following elements which can be modified.
#'Most users will not need to modify any of these, except possibly `max_iter`.
#'
#'kmax: Maximum number of factors. Defaults to twice the number of traits
#'
#'cond_num: Maximum allowable condition number for R. Defaults to 1e5.
#'
#'max_iter: Maximum number of iterations. Defaults to 1000.
#'
#'extrapolate: Passed to flashier, defaults to TRUE
#'
#'ebnm_fn_F and ebnm_fn_L: Prior families for factors and loadings. Defaults to point-normal.
#'
#'init_fn: Flashier initialization function.
#'
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
                    method = c("fixed_factors", "random_effect"),
                    mode = c("z-score", "b-std"),
                    no_wrapup = FALSE){


  ## process parameters
  default_params <- gfa_default_parameters()
  for(n in names(default_params)){
    if(is.null(params[[n]])) params[[n]] <- default_params[[n]]
  }
  for(n in names(params)){
    if(! n %in% names(default_params)) stop("Unknown parameter ", n, " provided.")
  }
  mode <- match.arg(mode)
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
  dat <- gfa_set_data(Y = Z_hat, scale = scale, R = R, params = params, mode = mode)

  method <- match.arg(method)

  if(is.null(dat$params$kmax)) dat$params$kmax <- 2*dat$p


  if(is.null(dat$R)){
    message("R is not supplied, fitting assuming full independence")
    fit <- fit_gfa_noR(dat)
    method <- "noR"
  }else if(method == "fixed_factors"){
    fit <- fit_gfa_ff(dat)
  }else if(method == "random_effect"){
    fit <- fit_gfa_re(dat)
  }

  #fit$method <- method
  ## wrap up
  if(is.null(fit$flash_fit$maxiter.reached) & !no_wrapup){
    fit <- fit %>% flash_nullcheck(remove = TRUE, tol = -Inf) # this will only remove 0 fators
    fit <- gfa_duplicate_check(fit,
                               dim = 2,
                               check_thresh = params$duplicate_check_thresh)
    ret <- gfa_wrapup(fit, method = method,
                      scale = dat$scale, nullcheck = TRUE)
    ret$mode <- mode
    ret$R <- dat$R
    ret$params <- dat$params
  }else{
    ret <- list(fit = fit,
                params = dat$params, scale = dat$scale,
                mode = mode, R = dat$R)
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
  #nf <- nmax # Number of fixed factors
  #Fitting
  # randomly initialize A
  fixed_ebnm = flash_ebnm(prior_family = "normal",
                          g_init = ashr::normalmix(pi = c(1),mean = c(0),sd = c(1)),
                          fix_g = TRUE)
  A_rand <- matrix(rnorm(n = dat$n * nmax, mean = 0, sd = 1), nrow = dat$n, ncol = nmax)

  #First initialize flash objects

  fit <-  flash_init(data = dat$Y, S = sqrt(lambda_min), var_type = 2) %>%
    flash_greedy(Kmax = dat$params$kmax,
                 init_fn = dat$params$init_fn,
                 ebnm_fn = list(dat$params$ebnm_fn_L, dat$params$ebnm_fn_F) )

  k <- fit$n_factors
  fit <- fit %>%
    flash_factors_init(., init = list(A_rand, W), ebnm_fn = fixed_ebnm) %>%
    flash_factors_fix(., kset = k + (1:nmax), which_dim = "factors") %>%
    flash_backfit(maxiter = dat$params$max_iter,
                  extrapolate=dat$params$extrapolate)
  fit$method <- "fixed_factors"
  return(fit)

}

fit_gfa_re <- function(dat){
  fit <-  flash_init(data = dat$Y, S = 0, var_type = 2, re_cov = dat$R)
  # fit <- flashier:::update_random_effect(fit$flash_fit) %>%
  #   flashier:::init.tau()
  # fit <- flashier:::set.obj(fit, flashier:::calc.obj(fit))
  fit <- fit %>% flash_greedy(Kmax = dat$params$kmax,
                              init_fn = dat$params$init_fn,
                              ebnm_fn = list(dat$params$ebnm_fn_L, dat$params$ebnm_fn_F) )

  fit <- fit %>% flash_backfit(maxiter = dat$params$max_iter,
                               extrapolate=dat$params$extrapolate,
                               verbose = dat$params$verbose)
  fit$method <- "random_effect"
  return(fit)
}

fit_gfa_noR <- function(dat){
  #First initialize flash objects
  fit <-  flash_init(data = dat$Y, S = dat$S, var_type = 2)

  # Add factors
  fit <- fit %>%
    flash_greedy(Kmax = dat$params$kmax,
                 init_fn = dat$params$init_fn,
                 ebnm_fn = list(dat$params$ebnm_fn_L, dat$params$ebnm_fn_F) ) %>%
    flash_backfit(maxiter = dat$params$max_iter, extrapolate = dat$params$extrapolate)
  fit$method <- "noR"
  return(fit)
}
