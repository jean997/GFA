#'@export
disc <- function(F_hat, F_true, lambda=0.9){
  b_est <- apply(cor(F_hat, F_true), 1, function(x){which.max(abs(x))})
  b_true <- apply(cor(F_hat, F_true), 2, function(x){which.max(abs(x))})
  b_est_v <-apply(cor(F_hat, F_true), 1, function(x){max(abs(x))})
  b_true_v <- apply(cor(F_hat, F_true), 2, function(x){max(abs(x))})

  matches <- which(b_true[b_est] == seq(length(b_est)))
  M <- data.frame(est_ix = matches, true_ix = b_est[matches], value = b_est_v[matches]) %>%
    filter(value >= lambda)
  unmatched_est_factors <- which(!seq_along(b_est) %in% M$est_ix)
  unmatched_true_factors <- which(!seq_along(b_true) %in% M$true_ix)
  R <- list("M"=M, u_est = unmatched_est_factors, u_true = unmatched_true_factors)
  return(R)
}
