
#'@export
est_R_pairwise <- function(B_hat, S_hat, N, p_val_thresh = 0.2, z_scores=FALSE){
  ntrait <- ncol(B_hat)
  nvar <- nrow(B_hat)
  stopifnot(nrow(S_hat) == nvar & ncol(S_hat)==ntrait)
  stopifnot(length(N)==ntrait)

  keep <- 2*pnorm(-abs(B_hat/S_hat)) > p_val_thresh
  if(z_scores){
    A <- B_hat/S_hat
  }else{
    A <- t( (1/sqrt(N)) *t(B_hat/S_hat))
  }

  R_df <- expand.grid(T1 = seq(ntrait), T2 = seq(ntrait)) %>%
        filter(T1 <= T2)
  r <- map2_dbl(R_df$T1, R_df$T2, function(i, j){
      if(i == j) return(1)
      ix <- which(keep[,i]==TRUE & keep[,j]==TRUE)
      cor(A[ix,i], A[ix,j])
  } )
  R_df$r <-r
  R <- reshape2::acast(R_df, T1~T2, value.var = "r")
  R_lower <- reshape2::acast(R_df, T2~T1, value.var = "r")
  R[lower.tri(R)] <- R_lower[lower.tri(R)]
  eig_R <- eigen(R)
  if(any(eig_R$values < 0)){
    eig_R$values <- pmin(eig_R$values, 0)
    R <- eiag_R$vectors %*% diag(eig_R$values) %*% t(eig_R$vectors)
  }
  return(R)
}
