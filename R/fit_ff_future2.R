
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
fit_ff_update <- function(Z_hat, B_std, N, R, kmax, ridge_penalty = 0,
                   min_ev = 1e-3, max_lr_percent = 1, lr_zero_thresh = 1e-10,
                   max_iter = 1000,
                   extrapolate = TRUE,
                   ebnm_fn_F = as.ebnm.fn(prior_family = "point_normal", optmethod = "nlm"),
                   ebnm_fn_L = as.ebnm.fn(prior_family = "point_normal", optmethod = "nlm"),
                   init_fn = flashier::init.fn.default,
                   fixed_truncate = Inf, duplicate_check_thresh = 0.5){


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
  if(missing(kmax)) kmax <- 2*ntrait


  if(missing(R)){
    warning("R is not supplied, fitting assuming full independence")
    if(mode == "zscore"){
      S <- 1
    }else if(mode == "std" ){
      S = t( (1/sqrt(N)) * t(matrix(1, nrow=nvar, ncol=ntrait)))
    }

    #First initialize flash objects
    fit <-  flash.init(data = Y, S = S, var.type = 2)

    # Add factors
    fit <- fit %>%
      flash.add.greedy(Kmax = kmax, init.fn = init_fn, ebnm.fn = list(ebnm_fn_L, ebnm_fn_F) ) %>%
      flash.backfit(maxiter = max_iter, extrapolate = extrapolate)
    if(is.null(fit$flash.fit$maxiter.reached)){
      fit <- gfa_duplicate_check(fit, dim = 2, check_thresh = duplicate_check_thresh)
      ret <- gfa_wrapup2(fit, nullcheck = TRUE)
    }else{
      ret <- list(fit = fit, fixed_ix = c())
    }
    return(ret)
  }

  stopifnot(nrow(R) == ntrait & ncol(R) == ntrait)
  if(mode == "std"){
    stopifnot(length(N) == ntrait)
    R <- diag(1/sqrt(N)) %*% R %*% diag(1/sqrt(N))
  }
  eS <- eigen(R)
  if(!all(eS$values >  min_ev)) stop("All eigenvalues of R must be greater than", min_ev)
  if(all(eS$values == 1)){
    warning("R appears to be the identity, you should rerun this command without it.")
    fit <- NULL
    Y <- NULL
    F_hat <- NULL
    L_hat <- NULL
    ret <- list(fit=fit, Y = Y, L_hat = L_hat, F_hat = F_hat)
    return(ret)
  }

  eS$values <- (eS$values + ridge_penalty)/(1 + ridge_penalty)

  lambda_min <- eS$values[ntrait]
  vals <- eS$values - lambda_min
  if(max_lr_percent < 1){
    vv <- cumsum(vals)/sum(vals)
    nmax <- min(which(vv > max_lr_percent)) - 1
  }else{
    nmax <- ntrait - 1
  }
  ix <- which(vals[seq(nmax)] > lr_zero_thresh)
  if(length(ix) == 0) stop("Is supplied R diagonal or close to diagonal?")

  W <- eS$vectors[,ix, drop = FALSE] %*% diag(sqrt(vals[ix]), ncol = length(ix))
  nf <- length(ix) # Number of fixed factors

  #Fitting
  # randomly initialize A
  if(fixed_truncate == Inf){
    fixed_ebnm = as.ebnm.fn(prior_family = "normal",
                            g_init = ashr::normalmix(pi = c(1), mean = c(0), sd = c(1)),
                            fix_g = TRUE)
    A_rand <- matrix(rnorm(n=nvar*nf, mean = 0, sd = 1), nrow=nvar, ncol=nf)
  }else{
    stop("Not implemented yet.")
    fixed_pfam = prior.unimodal.symmetric(scale = 1, g_init=fixed_g, fix_g = TRUE)
    A_rand <- matrix(runimix(n=nvar*nf, g = fixed_g), nrow=nvar, ncol=nf)
    A2_rand <- matrix(rnorm(n=nvar*nf), nrow=nvar, ncol=nf)
  }
  #First initialize flash objects

  fit <-  flash.init(data = Y, S = sqrt(lambda_min), var.type = 2) %>%
    flash.add.greedy(Kmax = kmax, init.fn = init_fn, ebnm.fn = list(ebnm_fn_L, ebnm_fn_F) )

  #Next add in fixed factors.
  n <- fit$n.factors
  fit <- fit %>%
         flash.init.factors(., init = list(A_rand, W), ebnm.fn = fixed_ebnm) %>%
         flash.fix.factors(., kset = n + (1:nf), mode=2) %>%
         flash.backfit(maxiter = max_iter, extrapolate=extrapolate)
  if(is.null(fit$flash.fit$maxiter.reached)){
    fit <- flash.nullcheck(remove = TRUE)
    fit <- gfa_duplicate_check(fit, dim = 2, check_thresh = duplicate_check_thresh)
    ret <- gfa_wrapup2(fit, nullcheck = TRUE)
  }else{
    ret <- list(fit = fit)
  }
  return(ret)
}

#'@export
gfa_wrapup2 <- function(fit, nullcheck = TRUE){
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
gfa_rebackfit2 <- function(fit, extrapolate = FALSE, maxiter){
  fit <- fit %>% flash.backfit(maxiter = maxiter, extrapolate = extrapolate)
  if(is.null(fit$flash.fit$maxiter.reached)){
    ret <- gfa_wrapup2(fit, nullcheck = TRUE)
  }else{
    ret <- list(fit = fit)
  }
  return(ret)
}

#'@export
est_L_flash2 <- function(Z_hat, fit, tol = 1e-5){
  flash_fit <- fit$flash.fit
  n_new <- nrow(Z_hat)
  s <- 1/sqrt(ff.tau(flash_fit))

  nfct <- fit$n.factors
  Lrand <- matrix(rnorm(n = n_new*nfct), nrow = n_new)
  g_ebnm <- list()
  fit_new <- flash.init(data = Z_hat, S = s, var.type = NULL)
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
}



