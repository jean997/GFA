
#'@title Fit with fixed factors
#'@param B_hat Matrix of (non-standardized) effect estimates SNPs by traits
#'@param S_hat Matrix of standard errors SNPs by traits
#'@param N Vector of sample sizes length equal to number of traits
#'@param R Estimated residual correlation matrix
#'@param kmax Maximum number of factors
#'@param zero_thresh Threshold for setting eigenvalues of R to zero
#'@return A list with elements fit, B_hat, L_hat, F_hat
#'@export
fit_ff <- function(B_hat, S_hat, N, R, kmax=100, zero_thresh = 1e-15, adjust=TRUE){

  n_var <- nrow(B_hat)
  n_trait <- ncol(B_hat)
  stopifnot(nrow(S_hat) == n_var & ncol(S_hat) == n_trait)
  #if(!adjust) stop("fit_ff only works with adjust right now.")
  if(!missing(N)) stopifnot(length(N) == n_trait)
  if(adjust & missing(N)) stop("To adjust please supply N.")

  R_eig <- eigen(R)
  R_eig$values[abs(R_eig$values) < zero_thresh] <- 0
  if(any(R_eig$values < 0))stop("R is not psd")

  Sigma <- diag(1/sqrt(N)) %*% R %*% diag(1/sqrt(N))
  lambda_min <- eigen(Sigma) %>%
    with(., min(values))
  Sig_new <- Sigma - diag(rep(lambda_min, n_trait))
  if(adjust){
    B_tilde = t( (1/sqrt(N)) *t(B_hat/S_hat))
    S_tilde = t( (1/sqrt(N)) * t(matrix(1, nrow=n_var, ncol=n_trait)))
  }else{
    B_tilde = B_hat
    S_tilde = S_hat
  }

  if(all(Sig_new ==0)){
    fit <- NULL
    B_hat <- NULL
    F_hat <- NULL
    L_hat <- NULL
  }else{
    Sig_eig <- eigen(Sig_new)
    V <- Sig_eig$vectors[, -n_trait]
    W <- V %*% diag(sqrt(Sig_eig$values[-n_trait]))

    # randomly initialize A
    A_rand <- matrix(rnorm(n=n_var*(n_trait-1)), nrow=n_var, ncol=(n_trait-1))

    #First add some greedy factors but don't backfit
    fit <-  flash.init(B_tilde, S = sqrt(lambda_min), var.type = 2) %>%
      flash.add.greedy(Kmax = n_trait, init.fn = init.fn.default )
    #Next add in fixed factors. Use sequential mode for backfit
    n <- fit$n.factors
    fit <- fit %>%
      flash.init.factors(., EF = list(A_rand, W), prior.family = prior.normal(scale= 1)) %>%
      flash.fix.loadings(., kset = n + 1:(n_trait-1), mode=2) %>%
      flash.backfit(method = "sequential")

    F_hat <- fit$loadings.pm[[2]][,1:n]
    L_hat <- fit$loadings.pm[[1]][, 1:n]
    fixed_ix <- n + (1:(n_trait-1))
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
