
#'@title Genetic Factor Analysis
#'@param Z_hat A matrix of z-scores with rows for variants and columns for traits.
#'@param N Vector of sample sizes length equal to number of traits. If not provided,
#'N will default to the vector of 1's and the factor matrix will be returned on the "z-score scale".
#'@param B_hat A matrix of GWAS effect estimates. B_hat is an alternative to Z_hat (only provide one of these).
#'If using B_hat, method must be "random_effect".
#'@param S If using B_hat, provide the corresponding matrix of standard errors.
#'@param R Estimated residual correlation matrix. This can be produced for example using R_ldsc, R_pt, or R_ldsc_quick.
#'@param params List of parameters. For most users this can be left at the default values. See details.
#'@return A list with elements fit, Y, L_hat, F_hat
#'@details
#'
#'The params list includes:
#'
#'kmax, the maximum number of factors
#'
#'max_iter: maximum number of iterations
#'
#'min_ev, max_lr_percent, lr_zero_thresh: We are approximating the matrix R as VDV^T + lambda I
#'where lambda is the smallest eigenvalue of R. min_ev specifies the smallest acceptable value of lambda.
#'max_lr_percent specifies the proportion of variance of (R - lambda I) that will be retained. This defaults
#'to 1 making the approximation exact. Eigenvalues of (R-lambda I) that are less than lr_zero_thresh will
#'be set to zero.
#'
#'extrapolate: value of extrapolate passed to flash_backfit
#'
#'ebnm_fn_F, ebnm_fn_L
#'
#'init_fn
#'
#'
#'duplicate_check_thresh
#'@export
gfa_fit <- function(Z_hat = NULL,
                    N = NULL,
                    B_hat = NULL,
                    S = NULL,
                    R = NULL,
                    params = gfa_default_parameters(),
                    method = c("fixed_factors", "random_effect"),
                    no_wrapup = FALSE, override = FALSE){


  ## process parameters
  default_params <- gfa_default_parameters()
  for(n in names(default_params)){
    if(is.null(params[[n]])) params[[n]] <- default_params[[n]]
  }
  for(n in names(params)){
    if(! n %in% names(default_params)) stop("Unknown parameter ", n, " provided.")
  }

  ## process inputs
  check_args_znbs(is.null(Z_hat), is.null(N), is.null(B_hat), is.null(S))
  if(is.null(Z_hat)){
    dat <- gfa_set_data(Y = B_hat, scale = NULL, S = S, R = R, params = params)
  }else{
    dat <- gfa_set_data(Y = Z_hat, scale = sqrt(N), S = NULL, R = R, params = params)
  }

  method <- match.arg(method)

  if(is.null(params$kmax)) params$kmax <- 2*dat$p


  if(is.null(dat$R)){
    message("R is not supplied, fitting assuming full independence")
    fit <- fit_gfa_noR(dat)
  }else if(method == "fixed_factors"){
    fit <- fit_gfa_ff(dat)
  }else if(method == "random_effect"){
    fit <- fit_gfa_re(dat)
  }

  ## wrap up
  if(is.null(fit$flash_fit$maxiter.reached) & !no_wrapup){
    fit <- fit %>% flash_nullcheck(remove = TRUE)
    fit <- gfa_duplicate_check(fit, dim = 2, check_thresh = params$duplicate_check_thresh)
    ret <- gfa_wrapup(fit, scale = dat$scale, nullcheck = FALSE)
    ret$params <- dat$params
    ret$scale <- dat$scale
  }else{
    ret <- list(fit = fit, params = dat$params, scale = dat$scale)
  }
  return(ret)
}

fit_gfa_ff <- function(dat){
  eS <- eigen(dat$R)
  lambda_min <- eS$values[dat$p]
  vals <- eS$values - lambda_min
  if (dat$params$max_lr_percent < 1) {
    vv <- cumsum(vals)/sum(vals)
    nmax <- min(which(vv > params$max_lr_percent))
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
  A_rand <- matrix(rnorm(n = nvar * nf, mean = 0, sd = 1), nrow = nvar, ncol = nf)

  #First initialize flash objects

  fit <-  flash_init(data = dat$Y, S = sqrt(lambda_min), var_type = 2) %>%
    flash_greedy(Kmax = dat$params$kmax,
                 init_fn = dat$params$init_fn,
                 ebnm_fn = list(dat$params$ebnm_fn_L, dat$params$ebnm_fn_F) )
  return(fit)

}

fit_gfa_re <- function(dat){
  fit <-  flash_init(data = dat$Y, S = 0, var_type = 2, re_cov = dat$R)
  fit <- flashier:::update_random_effect(fit$flash_fit) %>%
    flashier:::init.tau()
  fit <- flashier:::set.obj(fit, flashier:::calc.obj(fit))
  fit <- fit %>% flash_greedy(Kmax = dat$params$kmax,
                              init_fn = dat$params$init_fn,
                              ebnm_fn = list(dat$params$ebnm_fn_L, dat$params$ebnm_fn_F) )

  fit <- fit %>% flash_backfit(maxiter = dat$params$max_iter,
                               extrapolate=dat$params$extrapolate,
                               verbose = dat$params$verbose)
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
}
