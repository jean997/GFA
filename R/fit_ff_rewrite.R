
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
#'@param S_NULL If TRUE, will fit EBMF without supplying S.
#'@return A list with elements fit, B_hat, L_hat, F_hat
#'@export
fit_ff_new <- function(Z_hat, B_std, N, R, kmax,
                   zero_thresh = 1e-15, method="sequential",
                   max_ev_percent = 1, S_NULL = FALSE, reverse_fit_order = FALSE){

  if(!missing(Z_hat) & !missing(B_std)) stop("Please supply only one of Z_hat and B_std")
  if(!missing(B_std) & missing(N)) stop("If using B_std, N is required.")

  if(!missing(Z_hat)){
    mode <- "zscore"
    B_hat <- Z_hat
  }else{
    mode <- "std"
    B_hat <- B_std
  }
  ntrait <- ncol(B_hat)
  nvar <- nrow(B_hat)

  if(missing(R)){
    if(missing(kmax)) kmax <- 2*ntrait
    warning("R is not supplied, fitting assuming full independence")
    if(mode == "zscore" & !S_NULL){
      fit <- flash.init(data=B_hat, S = 1,  var.type=2)
    }else if(mode == "std" & !S_NULL){
      S = t( (1/sqrt(N)) * t(matrix(1, nrow=nvar, ncol=ntrait)))
      fit <- flash.init(data=B_hat, S = S,  var.type=2)
    }else{
      fit <- flash.init(data=B_hat, var.type=2)
    }
    fit <- fit %>%
           flash.add.greedy(Kmax = kmax, init.fn = init.fn.softImpute,
                              prior.family = prior.point.normal()) %>%
           flash.backfit() %>%
           flash.nullcheck()

    if(fit$n.factors > 0){
      F_hat <- fit$loadings.pm[[2]]
      L_hat <- fit$loadings.pm[[1]]
      B_hat <- fitted(fit)
      c <- colSums(F_hat^2)
      if(any(c==0)){
        i <- which(c==0)
        F_hat <- F_hat[,-i]
        L_hat <- L_hat[,-i]
      }
    }else{
      B_hat <- F_hat  <- L_hat <- NULL;
    }
    ret <- list(fit=fit, B_hat = B_hat, L_hat = L_hat, F_hat = F_hat)
    return(ret)
  }

  stopifnot(nrow(R) == ntrait & ncol(R) == ntrait)
  if(mode == "std"){
    stopifnot(lenght(N) == ntrait)
    R <- diag(1/sqrt(N)) %*% R %*% diag(1/sqrt(N))
  }
  eS <- eigen(R)
  eS$values[abs(eS$values) < zero_thresh] <- 0
  if(any(eS$values < 0))stop("R is not psd")
  if(all(eS$values == 0)) stop("All eigenvectors equal to zero.")
  if(all(eS$values == 1)){
    warning("R appears to be the identity, you should rerun this command without it.")
    fit <- NULL
    B_hat <- NULL
    F_hat <- NULL
    L_hat <- NULL
    ret <- list(fit=fit, B_hat = B_hat, L_hat = L_hat, F_hat = F_hat)
    return(ret)
  }
  if(max_ev_percent < 1){
    vv <- cumsum(eS$values)/sum(eS$values)
    nmax <- which.min(vv < max_ev_percent)-1
  }else{
    nmax <- sum(eS$values > 0)
  }

  lambda_min <- eS$values[nmax]
  V <- eS$vectors[,1:(nmax-1)]
  W <- V %*% diag(sqrt(eS$values[1:(nmax-1)]-lambda_min))
  nf <- nmax-1
  if(missing(kmax)) kmax <- ntrait
  #Fitting
  # randomly initialize A
  A_rand <- matrix(rnorm(n=nvar*nf), nrow=nvar, ncol=nf)
  gg <- ashr::normalmix(pi = c(1), mean = c(0), sd = c(1))
  if(!reverse_fit_order){
    #First add some greedy factors but don't backfit
    if(S_NULL){
        fit <-  flash.init(B_hat, var.type = 2) %>%
          flash.add.greedy(Kmax = kmax, init.fn = init.fn.default )
    }else{
      fit <-  flash.init(B_hat, S = sqrt(lambda_min), var.type = 2) %>%
        flash.add.greedy(Kmax = kmax, init.fn = init.fn.default )
    }
    #Next add in fixed factors. Use sequential mode for backfit
    n <- fit$n.factors
    fit <- fit %>%
           flash.init.factors(., EF = list(A_rand, W),
                              prior.family = prior.normal(scale= 1,
                              g_init=gg, fix_g = TRUE)) %>%
           flash.fix.loadings(., kset = n + 1:nf, mode=2) %>%
           flash.backfit(method = method)
    F_hat <- fit$loadings.pm[[2]][,1:n, drop=FALSE]
    L_hat <- fit$loadings.pm[[1]][, 1:n, drop=FALSE]
    fixed_ix <- n + (1:nf)
  }else{
    if(S_NULL){
      fit <-  flash.init(B_hat, var.type = 2) %>%
        flash.add.greedy(Kmax = kmax, init.fn = init.fn.default )
    }else{
      fit <-  flash.init(B_hat, S = sqrt(lambda_min), var.type = 2) %>%
        flash.add.greedy(Kmax = kmax, init.fn = init.fn.default )
    }
    #first add in fixed factors. Use sequential mode for backfit
    fit <- fit %>%
      flash.init.factors(., EF = list(A_rand, W),
                         prior.family = prior.normal(scale= 1,
                                                     g_init=gg, fix_g = TRUE)) %>%
      flash.fix.loadings(., kset = 1:nf, mode=2) %>%
      flash.add.greedy(Kmax = kmax, init.fn = init.fn.default ) %>%
      flash.backfit(method = method)
    F_hat <- fit$loadings.pm[[2]][,-(1:nf), drop=FALSE]
    L_hat <- fit$loadings.pm[[1]][, -(1:nf), drop=FALSE]
    fixed_ix <- 1:nf
  }
  B_hat <- fitted(fit) -
    with(fit, loadings.pm[[1]][, fixed_ix]%*%diag(loadings.scale[fixed_ix])%*% t(loadings.pm[[2]][, fixed_ix]))
  c <- colSums(F_hat^2)
  if(any(c==0)){
      i <- which(c==0)
      F_hat <- F_hat[,-i]
      L_hat <- L_hat[,-i]
  }
  ret <- list(fit=fit, B_hat = B_hat, L_hat = L_hat, F_hat = F_hat)
  return(ret)
}
