#'@export
sparse_svd <- function(dat, z_score_thresh){
  Z <- with(dat, beta_hat/se_beta_hat)
  Z[abs(Z) < z_score_thresh] <- 0
  svd(Z)
}
