#' @title Minimum norm distance between true and estimated factors
#' @description Given two matrices of factors (true and estimated), find the rotation
#'  that minimizes the Frobenius norm of the difference between the true factors and
#'  rotated estimate.
#' @param f_true Matrix of true factors (M x K1)
#' @param f_hat Matrix of estimated factors (M x K2)
#' @param single_trait_thresh Threshold to identify single trait factors. A factor is
#' considered a single trait factor if the maximum absolute value of its entries
#' is greater than this threshold after normalizing the factor to have unit norm.
#' Single trait factors are removed from both f_true and f_hat before matching.
#' Default is 0.98.
#' @param return_Q Logical. If TRUE, return the optimal rotation matrix.
#' @return A list with the following elements:
#' \item{solution}{A data frame with the following columns:
#' \itemize{
#' \item true_ix: Index of the true factor
#' \item est_ix: Index of the matching estimated factor
#' \item max_true_val: Maximum absolute value of the true factor
#' \item max_hat_val: Maximum absolute value of the estimated factor
#' \item penalty: The squared Frobenius norm penalty for the matched pair
#' \item match_score: The matching score (absolute inner product) for the matched pair
#' }}
#' \item{frob_n}{The Frobenius norm of the difference between the best matching
#' factors.}
#' \item{Q}{(optional) The optimal matching matrix. Returned if return_Q = TRUE.}
#' \item{best_est}{(optional) The best matching estimated factors. Returned if return_Q = TRUE.}
#' \item{best_true}{(optional) The true factors corresponding to the best matching. Returned if return_Q = TRUE.}
#'@export
min_norm <- function(f_true, f_hat,
                     single_trait_thresh = 0.98,
                     return_Q = FALSE){

  M <- nrow(f_true)

  stopifnot(nrow(f_hat) == M)
  f_hat <- norm_cols(f_hat)$A
  f_hat_orig <- f_hat
  hat_ix <- seq(ncol(f_hat))

  f_true <- norm_cols(f_true)$A
  f_true_orig <- f_true
  true_ix <- seq(ncol(f_true))

  # remove single trait factors from f_true
  sol_single <- NULL
  true_single <- find_single_trait(f_true, single_trait_thresh)
  if(length(true_single) > 0){
    cat("There are ", length(true_single), " single trait factors in f_true\n")
    f_true <- f_true_orig[,-true_single, drop = FALSE]
    true_ix <- true_ix[-true_single]
    sol_single <- data.frame(true_ix =  true_single,
                             est_ix = NA,
                             val = NA)
  }

  # remove single trait factors from f_hat
  hat_single <- find_single_trait(f_hat, single_trait_thresh)
  if(length(hat_single) > 0){
    cat("There are ", length(hat_single), " single trait factors in f_hat\n")
    f_hat <- f_hat_orig[,-hat_single, drop = FALSE]
    hat_ix <- hat_ix[-hat_single]
    if(is.null(sol_single)){
      sol_single <- data.frame(true_ix =  NA,
                               est_ix = hat_single,
                               val = NA)
    }else{
      sol_single <- bind_rows(sol_single,
                              data.frame(true_ix =  NA,
                                         est_ix = hat_single,
                                         val = NA))
    }
  }


  ## no factors estimated but there are true factors
  if(length(hat_ix) == 0  & length(true_ix) > 0){
    solution <- data.frame(true_ix = true_ix,
                           est_ix = NA,
                           penalty = 1,
                           match_score = 0,
                           max_hat_val = NA,
                           max_true_val = NA)
    solution$max_true_val <- apply(abs(f_true), 2, max)
    if( !is.null(sol_single)){
      solution <- bind_rows(solution, sol_single)
    }
    frob_n <- sqrt(length(true_ix))
    ret <- list(solution = solution,
                frob_n = frob_n)
    return(ret)
  }

  ## No true factors
  if(length(true_ix) ==0){
    if(length(hat_ix) == 0){
      ret <- list(solution = NULL, frob_n = 0)
    }else{
      solution = data.frame(true_ix = NA,
                            est_ix = hat_ix,
                            max_true_val = NA,
                            max_hat_val = NA,
                            penalty = 1,
                            match_score = 0)
      solution$max_hat_val <- apply(abs(f_hat), 2, max)
      if( !is.null(sol_single)){
        solution <- bind_rows(solution, sol_single)
      }
      ret <- list(solution = solution,
                  frob_n = sqrt(length(hat_ix)))
    }
    return(ret)
  }

  k <- ncol(f_true) - ncol(f_hat)
  if(k > 0){ # more true factors than estimated
    f_hat <- cbind(f_hat, matrix(0, nrow = M, ncol = k))
    hat_extra <- seq(k) + length(hat_ix)
  }else if(k < 0){ # more estimated than true
    f_true <- cbind(f_true, matrix(0, nrow = M, ncol = -k))
    true_extra <- seq(-k) + length(true_ix)
  }

  d <- t(f_true) %*% f_hat

  b <- lp.assign(abs(d), direction = "max")

  solution <- data.frame(true_ix = seq(ncol(d)),
                         est_ix = apply(b$solution,1, function(x){which(x > 0)}),
                         val = apply(abs(d)*b$solution,1, max))

  sgn <- purrr::map2(solution$true_ix, solution$est_ix, function(i,j){sign(d[i,j])}) %>% unlist()
  sgn[sgn == 0] <- 1
  Q <- t(b$solution*sgn)

  ## find factors where a single trait factor from f_hat is matched to a zero factor.
  ## f_hat should not be penalized for these
  fhq <- f_hat %*% Q
  z <- (f_true-fhq)^2
  solution$penalty <- colSums(z)


  if(k > 0){
    solution$est_ix[solution$est_ix %in% hat_extra] <- NA
  }else if(k < 0){
    solution$val[solution$tue_ix %in% true_extra] <- NA
    solution$true_ix[solution$true_ix %in% true_extra] <- NA
  }
  solution$true_ix <- true_ix[solution$true_ix]
  solution$est_ix <- hat_ix[solution$est_ix]
  if(!is.null(sol_single)){
    solution <- bind_rows(solution, sol_single)
  }



  solution$max_hat_val <- apply(abs(f_hat_orig[, solution$est_ix, drop = FALSE]), 2, max)
  solution$max_true_val <- apply(abs(f_true_orig[, solution$true_ix, drop = FALSE]), 2, max)

  frob_n <- sqrt(sum(z))

  solution <- arrange(solution, -1*val) %>%
              select(true_ix,
                     est_ix,
                     max_true_val,
                     max_hat_val,
                     penalty,
                     val) %>%
              rename(match_score = val)


  ret <- list(solution = solution,
              frob_n = frob_n)
  if(return_Q){
    ret$Q <- Q
    ret$best_est <-  f_hat %*% Q
    ret$best_true <- f_true
  }

  return(ret)
}



norm_cols <- function(A){
  w <- colSums(A^2)
  return(list("A" = t(t(A)/sqrt(w)), "w" = w))
}

find_single_trait <- function(f, st_thresh = 0.95){
  f_norm <- norm_cols(f)$A
  col_max <- apply(abs(f_norm), 2, max)
  i <- which(col_max > st_thresh)
  return(i)
}
