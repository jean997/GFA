
#'@title Fit with fixed factors
#'@param B_hat Matrix of (non-standardized) effect estimates SNPs by traits
#'@param S_hat Matrix of standard errors SNPs by traits
#'@param N Vector of sample sizes length equal to number of traits
#'@param R Estimated residual correlation matrix
#'@param kmax Maximum number of factors
#'@param zero_thresh Threshold for setting eigenvalues of R to zero
#'@return A list with elements fit, B_hat, L_hat, F_hat
#'@export
fit_ff <- function(B_hat, S_hat, N, R, kmax=100, zero_thresh = 1e-15, adjust=TRUE,
                   method="sequential", max_ev_percent = 1, S_NULL = FALSE){

  n_var <- nrow(B_hat)
  n_trait <- ncol(B_hat)
  stopifnot(nrow(S_hat) == n_var & ncol(S_hat) == n_trait)
  #if(!adjust) stop("fit_ff only works with adjust right now.")
  if(!missing(N)) stopifnot(length(N) == n_trait)
  if(adjust & missing(N)) stop("To adjust please supply N.")

  R_eig <- eigen(R)
  R_eig$values[abs(R_eig$values) < zero_thresh] <- 0
  if(any(R_eig$values < 0))stop("R is not psd")

  if(!all(N==1)){
    Sigma <- diag(1/sqrt(N)) %*% R %*% diag(1/sqrt(N))
    eS <- eigen(Sigma)
  }else{
    Sigma <- R
    eS <- R_eig
  }

  if(max_ev_percent < 1){
    stopifnot(S_NULL)
    vv <- cumsum(eS$value)/sum(eS$values)
    nmax <- which.min(vv < max_ev_percent)-1
    lambda_min <- min(eS$values)
    V <- eS$vectors[,1:nmax]
    W <- V %*% diag(sqrt(eS$values[1:nmax]-lambda_min))
    #X <-  V %*% diag(eS$values[1:nmax]) %*% t(V)
    #s <- matrix(rep(sqrt(diag(Sigma-X)), n_var), nrow=n_var, byrow = TRUE)
    nf <- nmax
  }else{

    lambda_min <- min(eS$values)
    V <- eS$vectors[, -n_trait]
    W <- V %*% diag(sqrt(eS$values[-n_trait]-lambda_min))
    s <- sqrt(lambda_min)
    nf <- n_trait -1
  }


  if(adjust){
    B_tilde = t( (1/sqrt(N)) *t(B_hat/S_hat))
  }else{
    B_tilde = B_hat
  }

  if(all(eS$values==0)){
    fit <- NULL
    B_hat <- NULL
    F_hat <- NULL
    L_hat <- NULL
  }else{

    # randomly initialize A
    A_rand <- matrix(rnorm(n=n_var*nf), nrow=n_var, ncol=nf)

    gg <- ashr::normalmix(pi = c(1), mean = c(0), sd = c(1))
    #First add some greedy factors but don't backfit

    if(S_NULL){
      fit <-  flash.init(B_tilde, var.type = 2) %>%
        flash.add.greedy(Kmax = n_trait, init.fn = init.fn.default )
      #Next add in fixed factors. Use sequential mode for backfit
      n <- fit$n.factors
      fit <- fit %>%
        flash.init.factors(., EF = list(A_rand, W),
                           prior.family = prior.normal(scale= 1,  g_init=gg, fix_g = TRUE)) %>%
        flash.fix.loadings(., kset = n + 1:nf, mode=2) %>%
        flash.backfit(method = method)
    }else{
      fit <-  flash.init(B_tilde, S = s, var.type = 2) %>%
        flash.add.greedy(Kmax = n_trait, init.fn = init.fn.default )
      #Next add in fixed factors. Use sequential mode for backfit
      n <- fit$n.factors
      fit <- fit %>%
        flash.init.factors(., EF = list(A_rand, W),
                           prior.family = prior.normal(scale= 1,  g_init=gg, fix_g = TRUE)) %>%
        flash.fix.loadings(., kset = n + 1:nf, mode=2) %>%
        flash.backfit(method = method)
    }
    F_hat <- fit$loadings.pm[[2]][,1:n, drop=FALSE]
    L_hat <- fit$loadings.pm[[1]][, 1:n, drop=FALSE]
    fixed_ix <- n + (1:nf)
    B_hat <- fitted(fit) -
      with(fit, loadings.pm[[1]][, fixed_ix]%*%diag(loadings.scale[fixed_ix])%*% t(loadings.pm[[2]][, fixed_ix]))
    c <- colSums(F_hat^2)
    if(any(c==0)){
      i <- which(c==0)
      F_hat <- F_hat[,-i]
      L_hat <- L_hat[,-i]
    }
  }
  ret <- list(fit=fit, B_hat = B_hat, L_hat = L_hat, F_hat = F_hat)
  return(ret)
}
