
#'@title Fit with fixed factors. This version is from March 2021 and explicitly only runs on
#'either z-scores or standardized effects. This function replaces both fit_ff and fit_plain.
#'@param Z_hat A matrix of z-scores
#'@param B_std A matrix of standardized effects. Supply only one of B_std or Z_hat
#'@param N Vector of sample sizes length equal to number of traits. N is only
#'needed if using B_std.
#'@param R Estimated residual correlation matrix
#'@param kmax Maximum number of factors. Defaults to 2*ntraits for plain or ntraits
#'if using the fixed factor correction.
#'@param max_ev_percent. Only the eigenvectors of R explaining the top percentage of
#'variance are retained. Defaults to 1, retaining all eigenvectors.
#'@param zero_thresh Threshold for setting eigenvalues of R to zero
#'@param S_inf A vector supplying a series of variance inflation factors for fitting ebmf.
#'S_inf should be decreasing and the last element should be 1.
#'@return A list with elements fit, Y, L_hat, F_hat
#'@export
fit_ff_prefit <- function(Z_hat, B_std, N, R, kmax,
                   zero_thresh = 1e-15, method="extrapolate",
                   max_ev_percent = 1,
                   S_inf = c(10, 1),
                   num_prefits = 1,
                   max_prefit_iter = 500,
                   min_var_ratio = 1){

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

  #Check if S_inf is valid
  if(missing(S_inf)){
    if(num_prefits == 0){
      S_inf <- c(1)
    }else{
      if(mode == "std") t <- 1/N
        else t <- 1
      max_var <- max(apply(Y, 2, var)/t)
      if(max_var < min_var_ratio){
        S_inf <- c(1)
      }else{
        max_v_inf <- max_var/min_var_ratio
        v_inf <- exp(seq(log(max_v_inf), 0, length.out = num_prefits+1))
        S_inf <- sqrt(v_inf)
      }
    }
  }

  if(missing(R)){
    if(missing(kmax)) kmax <- 2*ntrait
    warning("R is not supplied, fitting assuming full independence")
    if(mode == "zscore"){
      S <- 1
    }else if(mode == "std" ){
      S = t( (1/sqrt(N)) * t(matrix(1, nrow=nvar, ncol=ntrait)))
    }
    fits <- lapply(S_inf, function(s){flash.init(data=Y, S = s*S,  var.type=2)})
    fits[[1]] <- fits[[1]] %>%
           flash.add.greedy(Kmax = kmax, init.fn = init.fn.softImpute,
                              prior.family = prior.point.normal())

    if(length(S_inf) > 1){
      for(i in 2:length(S_inf)){
        fits[[i-1]] <- fits[[i-1]] %>%
                      flash.backfit(maxiter = max_prefit_iter)
        fits[[i]] <- fits[[i]] %>%
                     flash.init.factors(EF = fits[[i-1]]$flash.fit$EF, EF2 = fits[[i-1]]$flash.fit$EF2) %>%
                     flash.backfit() %>%
                     flash.add.greedy(Kmax = kmax, init.fn = init.fn.softImpute,
                           prior.family = prior.point.normal())
      }
    }
    fit <- fits[[length(S_inf)]]
    fit <- fit %>%
           flash.backfit() %>%
           flash.nullcheck(remove = TRUE)
    if(fit$n.factors > 0){
      F_hat <- fit$loadings.pm[[2]]
      L_hat <- fit$loadings.pm[[1]]
      Yhat <- fitted(fit)
      c <- colSums(F_hat^2)
      if(any(c==0)){
        i <- which(c==0)
        F_hat <- F_hat[,-i]
        L_hat <- L_hat[,-i]
      }
    }else{
      Yhat <- F_hat  <- L_hat <- NULL;
    }
    ret <- list(fit=fit, B_hat = Yhat, L_hat = L_hat, F_hat = F_hat)
    return(ret)
  }

  stopifnot(nrow(R) == ntrait & ncol(R) == ntrait)
  if(mode == "std"){
    stopifnot(length(N) == ntrait)
    R <- diag(1/sqrt(N)) %*% R %*% diag(1/sqrt(N))
  }
  eS <- eigen(R)
  eS$values[abs(eS$values) < zero_thresh] <- 0
  if(any(eS$values < 0))stop("R is not psd")
  if(all(eS$values == 0)) stop("All eigenvectors equal to zero.")
  if(all(eS$values == 1)){
    warning("R appears to be the identity, you should rerun this command without it.")
    fit <- NULL
    Y <- NULL
    F_hat <- NULL
    L_hat <- NULL
    ret <- list(fit=fit, Y = Y, L_hat = L_hat, F_hat = F_hat)
    return(ret)
  }
  if(max_ev_percent < 1){
    vv <- cumsum(eS$values)/sum(eS$values)
    nmax <- which.min(vv < max_ev_percent)-1
  }else{
    nmax <- sum(eS$values > 0)
  }

  lambda_min <- eS$values[nmax]
  V <- eS$vectors[,1:(nmax-1), drop = FALSE]
  W <- V %*% diag(sqrt(eS$values[1:(nmax-1)]-lambda_min), ncol = (nmax -1))
  nf <- nmax-1
  if(missing(kmax)) kmax <- ntrait
  #Fitting
  # randomly initialize A
  A_rand <- matrix(rnorm(n=nvar*nf), nrow=nvar, ncol=nf)
  gg <- ashr::normalmix(pi = c(1), mean = c(0), sd = c(1))

  #First initialize flash objects
  fits <-  lapply(S_inf, function(s){flash.init(data = Y, S = s*sqrt(lambda_min), var.type = 2)})

  fits[[1]] <- fits[[1]] %>%
               flash.add.greedy(Kmax = kmax, init.fn = init.fn.default )
  #Next add in fixed factors.
  n <- fits[[1]]$n.factors
  fits[[1]] <- fits[[1]] %>%
         flash.init.factors(., EF = list(A_rand, W),
                            prior.family = prior.normal(scale= 1,
                            g_init=gg, fix_g = TRUE)) %>%
         flash.fix.loadings(., kset = n + (1:nf), mode=2)

  fixed_ix <- n + (1:nf)
  if(length(S_inf) > 1){
    for(i in 2:length(S_inf)){
      fits[[i-1]] <- fits[[i-1]] %>%
                     flash.backfit(method = method, maxiter = max_prefit_iter)
      fits[[i]] <- fits[[i]] %>%
                   flash.init.factors(EF = fits[[i-1]]$flash.fit$EF, EF2 = fits[[i-1]]$flash.fit$EF2) %>%
                   flash.fix.loadings(., kset = fixed_ix, mode=2)
      fits[[i]]$flash.fit$is.zero <- fits[[i-1]]$flash.fit$is.zero
      fits[[i]] <- fits[[i]] %>%
                   flash.backfit(method = method, maxiter = max_prefit_iter) %>%
                   flash.add.greedy(Kmax = kmax, init.fn = init.fn.softImpute,
                         prior.family = prior.point.normal())
    }
  }
  fit <- fits[[length(S_inf)]]
  fit <- fit %>%
         flash.backfit(method = method) %>%
         flash.nullcheck(remove = FALSE)
  F_hat <- fit$loadings.pm[[2]][,-fixed_ix, drop=FALSE]
  L_hat <- fit$loadings.pm[[1]][, -fixed_ix, drop=FALSE]

  Yhat <- fitted(fit) -
    with(fit, loadings.pm[[1]][, fixed_ix, drop =FALSE]%*%diag(loadings.scale[fixed_ix], ncol = length(fixed_ix))%*% t(loadings.pm[[2]][, fixed_ix, drop = FALSE]))
  c <- colSums(F_hat^2)
  if(any(c==0)){
      i <- which(c==0)
      F_hat <- F_hat[,-i]
      L_hat <- L_hat[,-i]
  }
  ret <- list(fit=fit, B_hat = Yhat, L_hat = L_hat, F_hat = F_hat)
  return(ret)
}
