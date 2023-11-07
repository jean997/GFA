
#'@title Fit with fixed factors. This version is from March 2021 and explicitly only runs on
#'either z-scores or standardized effects. This function replaces both fit_ff and fit_plain.
#'@param Z_hat A matrix of z-scores
#'@param B_std A matrix of standardized effects. Supply only one of B_std or Z_hat
#'@param N Vector of sample sizes length equal to number of traits. N is only
#'needed if using B_std.
#'@param R Estimated residual correlation matrix
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
#'ridge_penalty
#'
#'extrapolate: value of extrapolate passed to flash_backfit
#'
#'ebnm_fn_F, ebnm_fn_L
#'
#'init_fn
#'
#'fixed_truncate
#'
#'duplicate_check_thresh
#'@export
gfa_fit <- function(Z_hat, B_std, N, R, params = gfa_default_parameters(),
                    no_wrapup = FALSE){


  default_params <- gfa_default_parameters()
  for(n in names(default_params)){
    if(is.null(params[[n]])) params[[n]] <- default_params[[n]]
  }
  for(n in names(params)){
    if(! n %in% names(default_params)) stop("Unknown parameter ", n, " provided.")
  }

  if(!missing(Z_hat) & !missing(B_std)) stop("Please supply only one of Z_hat and B_std")
  if(!missing(B_std) & missing(N)) stop("If using B_std, N is required.")
  if(!missing(Z_hat)){
    mode <- "zscore"
    Y <- Z_hat
  }else{
    mode <- "std"
    Y <- B_std
  }
  ntrait <- ncol(Y)
  nvar <- nrow(Y)
  if(is.null(params$kmax)) params$kmax <- 2*ntrait


  if(missing(R)){
    warning("R is not supplied, fitting assuming full independence")
    if(mode == "zscore"){
      S <- 1
    }else if(mode == "std" ){
      S = t((1/sqrt(N)) * matrix(1, nrow=ntrait, ncol=nvar))
    }

    #First initialize flash objects
    fit <-  flash_init(data = Y, S = S, var_type = 2)

    # Add factors
    fit <- fit %>%
      flash_greedy(Kmax = params$kmax,
                  init_fn = params$init_fn,
                  ebnm_fn = list(params$ebnm_fn_L, params$ebnm_fn_F) ) %>%
      flash_backfit(maxiter = params$max_iter, extrapolate = params$extrapolate)
    if(is.null(fit$flash_fit$maxiter.reached)){
      fit <- fit %>% flash_nullcheck(remove = TRUE)
      fit <- gfa_duplicate_check(fit, dim = 2, check_thresh = params$duplicate_check_thresh)
      ret <- gfa_wrapup(fit, nullcheck = TRUE)
      ret$params <- params
    }else{
      ret <- list(fit = fit, fixed_ix = c(), params = params)
    }
    return(ret)
  }

  stopifnot(nrow(R) == ntrait & ncol(R) == ntrait)
  if(mode == "std"){
    stopifnot(length(N) == ntrait)
    R <- diag(1/sqrt(N)) %*% R %*% diag(1/sqrt(N))
  }
  eS <- eigen(R)

  vals <- eS$values

  if(min(vals)/max(vals) < params$min_ev){
    stop("Ratio of smallest eigenvalue to largest eigenvalue is smaller than ", params$min_ev)
  }

  v <- sum((vals - min(vals))^2)/sum(vals^2)
  if(v < params$lr_zero_thresh){
    warning("R appears to be the identity or very close to the identity, fitting model assuming R = I.")
    if(mode == "std") res <- gfa_fit(B_std = B_std, N = N, params = params)
    else res <- gfa_fit(Z_hat = Z_hat, params = params)
    return(res)
  }

  vals <- (vals + params$ridge_penalty)/(1 + params$ridge_penalty)



  lambda_min <- vals[ntrait]
  vals <- vals - lambda_min
  if (params$max_lr_percent < 1) {
    vv <- cumsum(vals)/sum(vals)
    nmax <- min(which(vv > params$max_lr_percent))
  }else{
    nmax <- ntrait - 1
  }


  W <-eS$vectors[, seq(nmax), drop = FALSE] %*% diag(sqrt(vals[seq(nmax)]), ncol = nmax)
  nf <- nmax # Number of fixed factors

  #Fitting
  # randomly initialize A
  fixed_ebnm = flash_ebnm(prior_family = "normal",
                            g_init = ashr::normalmix(pi = c(1),mean = c(0),sd = c(1)),
                            fix_g = TRUE)
  A_rand <- matrix(rnorm(n = nvar * nf, mean = 0, sd = 1), nrow = nvar, ncol = nf)

  #First initialize flash objects

  fit <-  flash_init(data = Y, S = sqrt(lambda_min), var_type = 2) %>%
    flash_greedy(Kmax = params$kmax,
                  init_fn = params$init_fn,
                  ebnm_fn = list(params$ebnm_fn_L, params$ebnm_fn_F) )

  #Next add in fixed factors.
  n <- fit$n_factors
  fit <- fit %>%
         flash_factors_init(., init = list(A_rand, W), ebnm_fn = fixed_ebnm) %>%
         flash_factors_fix(., kset = n + (1:nf), which_dim = "factors") %>%
         flash_backfit(maxiter = params$max_iter,
                       extrapolate=params$extrapolate)

  if(is.null(fit$flash_fit$maxiter.reached) & !no_wrapup){
    fit <- fit %>% flash_nullcheck(remove = TRUE)
    fit <- gfa_duplicate_check(fit, dim = 2, check_thresh = params$duplicate_check_thresh)
    ret <- gfa_wrapup(fit, nullcheck = TRUE)
    ret$params <- params
  }else{
    ret <- list(fit = fit, params = params)
  }
  return(ret)
}

#'@export
gfa_wrapup <- function(fit, nullcheck = TRUE){
  if(nullcheck){
    fit <- fit %>% flash_nullcheck(remove = TRUE)
  }
  if(length(fit$flash_fit$fix.dim) == 0){
    fixed_ix <- rep(FALSE, ncol(fit$F_pm))
  }else{
    fixed_ix <- fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)})
  }

  if(any(fixed_ix)){
    fixed_ix <- which(fixed_ix)
    F_hat <- fit$F_pm[,-fixed_ix, drop=FALSE]
    L_hat <- fit$L_pm[, -fixed_ix, drop=FALSE]
  }else{
    F_hat <- fit$F_pm
    L_hat <- fit$L_pm
  }
  Yhat <- L_hat %*% t(F_hat)
  Lx <- colSums(L_hat^2)
  Fx <- colSums(F_hat^2)
  F_hat_scaled <- t(t(F_hat)/sqrt(Fx))
  L_hat_scaled <- t(t(L_hat)/sqrt(Lx))
  scale <- sqrt(Fx)*sqrt(Lx)
  ret <- list(fit=fit, B_hat = Yhat, L_hat = L_hat, F_hat = F_hat,
              scale = scale, fixed_ix = fixed_ix,
              F_hat_scaled = F_hat_scaled,
              L_hat_scaled = L_hat_scaled)
  return(ret)
}

#'@export
gfa_rebackfit <- function(fit, params){
  fit <- fit %>% flash_backfit(maxiter = params$max_iter,
                               extrapolate = params$extrapolate)
  if(is.null(fit$flash_fit$maxiter.reached)){
    ret <- gfa_wrapup(fit, nullcheck = TRUE)
  }else{
    ret <- list(fit = fit)
  }
  ret$params <- params
  return(ret)
}




