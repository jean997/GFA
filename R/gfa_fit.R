
#'@title Fit with fixed factors. This version is from March 2021 and explicitly only runs on
#'either z-scores or standardized effects. This function replaces both fit_ff and fit_plain.
#'@param Z_hat A matrix of z-scores
#'@param B_std A matrix of standardized effects. Supply only one of B_std or Z_hat
#'@param N Vector of sample sizes length equal to number of traits. N is only
#'needed if using B_std.
#'@param R Estimated residual correlation matrix
#'@param kmax Maximum number of factors. Defaults to 2*ntraits for plain or ntraits
#'if using the fixed factor correction.
#'@param min_ev,max_lr_percent,lr_zero_thresh See details
#'@param S_inf A vector supplying a series of variance inflation factors for fitting ebmf.
#'S_inf should be decreasing and the last element should be 1.
#'@param fixed_g Supply the prior for loadings of fixed factors.
#'@return A list with elements fit, Y, L_hat, F_hat
#'@details
#'
#'min_ev, max_lr_percent, lr_zero_thresh: We are approximating the matrix R as VDV^T + lambda I
#'where lambda is the smallest eigenvalue of R. min_ev specifies the smallest acceptable value of lambda.
#'max_lr_percent specifies the proportion of variance of (R - lambda I) that will be retained. This defaults
#'to 1 making the approximation exact. Eigenvalues of (R-lambda I) that are less than lr_zero_thresh will
#'be set to zero.
#'
#'@export
gfa_fit <- function(Z_hat, B_std, N, R, params = gfa_default_parameters()){


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
    fit <-  flash.init(data = Y, S = S, var.type = 2)

    # Add factors
    fit <- fit %>%
      flash.add.greedy(Kmax = params$kmax,
                       init.fn = params$init_fn,
                       ebnm.fn = list(params$ebnm_fn_L, params$ebnm_fn_F) ) %>%
      flash.backfit(maxiter = params$max_iter, extrapolate = params$extrapolate)
    if(is.null(fit$flash.fit$maxiter.reached)){
      fit <- fit %>% flash.nullcheck(remove = TRUE)
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
    nmax <- min(which(vv > params$max_lr_percent)) - 1
  }else{
    nmax <- ntrait - 1
  }
  ix <- which(vals[seq(nmax)] > params$lr_zero_thresh)

  W <-eS$vectors[, ix, drop = FALSE] %*% diag(sqrt(vals[ix]), ncol = length(ix))
  nf <- length(ix) # Number of fixed factors

  #Fitting
  # randomly initialize A
  if (params$fixed_truncate == Inf) {
    fixed_ebnm = as.ebnm.fn(prior_family = "normal",
                            g_init = ashr::normalmix(pi = c(1),mean = c(0),sd = c(1)),
                            fix_g = TRUE)
    A_rand <- matrix(rnorm(n = nvar * nf, mean = 0, sd = 1), nrow = nvar, ncol = nf)
  }else{
    stop("Not implemented yet.")
  }

  #First initialize flash objects

  fit <-  flash.init(data = Y, S = sqrt(lambda_min), var.type = 2) %>%
    flash.add.greedy(Kmax = params$kmax,
                     init.fn = params$init_fn,
                     ebnm.fn = list(params$ebnm_fn_L, params$ebnm_fn_F) )

  #Next add in fixed factors.
  n <- fit$n.factors
  fit <- fit %>%
         flash.init.factors(., init = list(A_rand, W), ebnm.fn = fixed_ebnm) %>%
         flash.fix.factors(., kset = n + (1:nf), mode=2) %>%
         flash.backfit(maxiter = params$max_iter,
                       extrapolate=params$extrapolate)

  if(is.null(fit$flash.fit$maxiter.reached)){
    fit <- fit %>% flash.nullcheck(remove = TRUE)
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
    fit <- fit %>% flash.nullcheck(remove = TRUE)
  }
  fixed_ix <- fit$flash.fit$fix.dim %>% sapply(., function(x){!is.null(x)})
  if(any(fixed_ix)){
    fixed_ix <- which(fixed_ix)
    F_hat <- fit$F.pm[,-fixed_ix, drop=FALSE]
    L_hat <- fit$L.pm[, -fixed_ix, drop=FALSE]
  }else{
    F_hat <- fit$F.pm
    L_hat <- fit$L.pm
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
  fit <- fit %>% flash.backfit(maxiter = params$max_iter,
                               extrapolate = params$extrapolate)
  if(is.null(fit$flash.fit$maxiter.reached)){
    ret <- gfa_wrapup(fit, nullcheck = TRUE)
  }else{
    ret <- list(fit = fit)
  }
  ret$params <- params
  return(ret)
}




