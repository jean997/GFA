gfa_greedy_init_softImpute <- function(flash, eig_sigma, seed = 666, maxiter = 1000, ...) {
  set.seed(seed)

  if (flashier:::get.dim(flash) > 2)
    stop("softImpute cannot be used with tensors.")

  if (inherits(flashier:::get.Y(flash), "lowrank"))
    stop("softImpute cannot be used with low-rank matrix representations.")

  if (inherits(flashier:::get.Y(flash), "lrps"))
    stop("softImpute cannot be used with low-rank plus sparse representations.")

  p <- length(eig_sigma$values)
  sqrt_inv_sigma <- eig_sigma$vectors %*% diag(1/sqrt(eig_sigma$values), nrow = p)
  tsqrt_sigma <- eig_sigma$vectors %*% diag(sqrt(eig_sigma$values), nrow = p)

  X <- residuals(flash)
  X <-  X %*% sqrt_inv_sigma

  si.res <- softImpute(X, rank.max = 1, maxit = maxiter, ...)

  v <- tsqrt_sigma %*% matrix(sqrt(si.res$d)*si.res$v, ncol = 1)
  v <- as.numeric(v)

  EF <- list(si.res$u * sqrt(si.res$d),v)

  return(EF)
}
