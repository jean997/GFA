#'@export
min_norm <- function(f_true, f_hat, l_true, l_hat,
                     single_trait_thresh = 0.95, thresh = 0,
                     return_Q = FALSE){

  M <- nrow(f_true)
  stopifnot(nrow(f_hat) == M)

  #pre-processing
  f_true <- norm_cols(f_true)$A
  f_hat <- norm_cols(f_hat)$A
  n_t <- ncol(f_true)
  n_h <- ncol(f_hat)

  col_max_true <- apply(abs(f_true), 2, max)
  col_max_hat <- apply(abs(f_hat), 2, max)
  hat_single <- true_single <- c()
  true_ix <- seq(n_t)
  hat_ix <- seq(n_h)

  if(any(col_max_true > single_trait_thresh)){
    i <- which(col_max_true > single_trait_thresh)
    true_single <- i
    cat("Removing ", length(i), " single trait factors from f_true\n")
    f_true <- f_true[,-i, drop = FALSE]
    true_ix <- true_ix[-i]
  }
  if(n_h == 0){
    solution <- data.frame(true_ix = true_ix, hat_ix = NA, val = 0 )
    if( length(true_single) > 0){
      single_df <- data.frame(true_ix =  true_single,
                              est_ix = NA,
                              val = NA)
      solution <- bind_rows(solution, single_df)
    }
    frob_n <- opt_frob_n <- length(true_ix)
    hat_ix <- c()
    ret <- list(solution = solution,
                frob_n = frob_n,
                true_ix = true_ix,
                hat_ix = hat_ix,
                opt_frob_n = opt_frob_n)
    if(!missing(l_hat)){
      ret$frob_n_l <- frob_n
    }
    return(ret)
  }

  if(any(col_max_hat > single_trait_thresh)){
    i <- which(col_max_hat > single_trait_thresh)
    hat_single <- i
    cat("Removing ", length(i), " single trait factors from f_hat\n")
    f_hat <- f_hat[,-i, drop = FALSE]
    hat_ix <- hat_ix[-i]
  }

  if(!missing(l_true) & !missing(l_hat)){
    stopifnot(ncol(l_true) == n_t & ncol(l_hat) == n_h)
    N <- nrow(l_true)
    stopifnot(nrow(l_hat) == N)
    l_true <- norm_cols(l_true)$A
    l_hat <- norm_cols(l_hat)$A
    l_true <- l_true[,true_ix]
    l_hat <- l_hat[,hat_ix]
    l <- TRUE
  }else{
    l <- FALSE
  }


  k <- ncol(f_true) - ncol(f_hat)
  if(k > 0){
    f_hat <- cbind(f_hat, matrix(0, nrow = M, ncol = k))
    hat_ix <- c(hat_ix, seq(k) + n_h)
    if(l){
      l_hat <- cbind(l_hat, matrix(0, nrow = N, ncol = k))
    }
  }else if(k < 0){
    f_true <- cbind(f_true, matrix(0, nrow = M, ncol = -k))
    true_ix <- c(true_ix, seq(-k) + n_t)
    if(l){
      l_true <- cbind(l_true, matrix(0, nrow = N, ncol = -k))
    }
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
  if(l){
    frob_n_l <- sum((l_true - l_hat%*%Q)^2)
  }
  if(k >= 0){
    opt_frob_n <- frob_n
  }else{
    # If more estimated factors than true, choose the optimal number
    opt_k <- n_t - length(true_single)
    opt_frob_n <- sum((f_true[, 1:opt_k] - (f_hat %*% Q)[, 1:opt_k])^2)
  }

  n <- ncol(d)


  solution <- solution %>% arrange(-val)
  solution$est_ix <- hat_ix[solution$est_ix]
  solution$true_ix <- true_ix[solution$true_ix ]
  solution$est_ix[solution$est_ix > n_h] <- NA
  solution$true_ix[solution$true_ix > n_t] <- NA

  if(length(hat_single) + length(true_single) > 0){
    single_df <- data.frame(true_ix = c(rep(NA, length(hat_single)), true_single),
                            est_ix = c(hat_single, rep(NA, length(true_single))),
                            val = NA)
    solution <- bind_rows(solution, single_df)
  }


  ret <- list(solution = solution,
              frob_n = frob_n,
              true_ix = true_ix,
              opt_frob_n = opt_frob_n,
              hat_ix = hat_ix)
  if(return_Q){
    ret$Q <- Q
  }
  if(l){
    ret$frob_n_l <- frob_n_l
  }
  return(ret)
}



norm_cols <- function(A){
  w <- colSums(A^2)
  return(list("A" = t(t(A)/sqrt(w)), "w" = w))
}
