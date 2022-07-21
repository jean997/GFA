#'@export
min_norm <- function(f_true, f_hat, single_trait_thresh = 0.95, thresh = 0){

  M <- nrow(f_true)
  stopifnot(nrow(f_hat) == M)

  #pre-processing
  f_true <- norm_cols(f_true)$A
  f_hat <- norm_cols(f_hat)$A
  col_max_true <- apply(abs(f_true), 2, max)
  col_max_hat <- apply(abs(f_hat), 2, max)
  hat_single <- true_single <- c()
  if(any(col_max_true > single_trait_thresh)){
    i <- which(col_max_true > single_trait_thresh)
    true_single <- i
    cat("Removing ", length(i), " single trait factors from f_true\n")
    f_true <- f_true[,-i]
  }
  if(any(col_max_hat > single_trait_thresh)){
    i <- which(col_max_hat > single_trait_thresh)
    hat_single <- i
    cat("Removing ", length(i), " single trait factors from f_hat\n")
    f_hat <- f_hat[,-i]
  }
  k <- ncol(f_true) - ncol(f_hat)
  if(k > 0){
    f_hat <- cbind(f_hat, matrix(0, nrow = M, ncol = k))
  }else if(k < 0){
    f_true <- cbind(f_true, matrix(0, nrow = M, ncol = -k))
  }




  d <- t(f_true) %*% f_hat
  d[abs(d) < thresh ] <- 0

  b <- lp.assign(abs(d), direction = "max")

  solution <- data.frame(true_ix = seq(ncol(d)),
                         est_ix = apply(b$solution,1, function(x){which(x > 0)}),
                         val = apply(abs(d)*b$solution,1, max))

  sgn <- purrr::map2(solution$true_ix, solution$est_ix, function(i,j){sign(d[i,j])}) %>% unlist()
  sgn[sgn == 0] <- 1
  Q <- t(b$solution*sgn)

  frob_n <- sum((f_true - f_hat%*%Q)^2)

  n <- ncol(d)
  if(k > 0){
    solution$est_ix[solution$est_ix > (n -k)] <- NA
  }else if(k  < 0){
    solution$true_ix[solution$true_ix > (n +k)] <- NA
  }

  solution <- solution %>% arrange(-val)
  solution$est_ix <- solution$est_ix + sapply(solution$est_ix, function(j){sum(hat_single <= j)})
  solution$true_ix <- solution$true_ix + sapply(solution$true_ix, function(j){sum(true_single <= j)})

  if(length(hat_single) + length(true_single) > 0){
    single_df <- data.frame(true_ix = c(rep(NA, length(hat_single)), true_single),
                            est_ix = c(hat_single, rep(NA, length(true_single))),
                            val = NA)
    solution <- bind_rows(solution, single_df)
  }


  return(list(solution = solution, val = frob_n))


}


norm_cols <- function(A){
  w <- colSums(A^2)
  return(list("A" = t(t(A)/sqrt(w)), "w" = w))
}
