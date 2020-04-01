
#'@export
sim_bh_simple <- function(rloadings, rfactors, rerror, S, k){
  nvar <- nrow(S)
  ntrait <- ncol(S)
  true_L <- replicate(n=k, rloadings(nvar))
  true_F <- replicate(n=k, rfactors(ntrait))
  true_E <- replicate(n=ntrait, rerror(nvar))
  beta <- true_L %*% t(true_F) + true_E
  beta_hat <- map_dfc(seq(ntrait),
                      function(i){rnorm(n=nvar, mean = beta[,i], sd = S[,i])}
  )
  mats <- list(beta_hat = as.matrix(beta_hat), se_hat = S, beta = beta,
               true_L = true_L, true_F = true_F, true_E = true_E)
  return(mats)
}

#'@title Simulate data from sparse factor model
#'@description Simulate B_hat from model B_hat = L^TF + Theta + E where elements of
#'E have E_ij ~ N(0, s_{ij}^2) and Cov(E_{i,.}) = SRS where R is M times M matrix of
#'sample correlation.
#'@param true_L Loadings matrix
#'@param true_F Factors matrix
#'@param true_Ttheta Theta
#'@param S standard errors
#'@param R Correlation matrix for rows of E
#'@export
sim_bh2<- function(true_L, true_F, true_Theta, S, R){
  M <- nrow(true_F)
  p <- nrow(true_L)
  K <- ncol(true_L)
  stopifnot(ncol(true_F) == K)
  stopifnot(nrow(true_Theta)==p & ncol(true_Theta) == M)
  stopifnot(nrow(S)==p & ncol(S) == M)
  stopifnot(nrow(R)==M & ncol(R) == M)

  true_E <- apply(S, 1, function(si){
      sigma <- diag(si) %*% R %*% diag(si)
      MASS::mvrnorm(n=1, mu = rep(0, M), Sigma = sigma)
  }) %>% t()
  beta_hat <- true_L %*% t(true_F) + true_Theta + true_E
  mats <- list(beta_hat = as.matrix(beta_hat), se_hat = S,
               true_L = true_L, true_F = true_F, true_E = true_E, true_Theta = true_Theta)
  return(mats)
}

#'@title Simulate data from sparse factor model
#'@description Simulate B_hat from model B_hat = L^TF + Theta + E where elements of
#'E have E_ij ~ N(0, s_{ij}^2) and Cov(E_{i,.}) = SRS where R is M times M matrix of
#'sample correlation.
#'@param true_L Loadings matrix
#'@param true_F Factors matrix
#'@param true_Ttheta Theta
#'@param S standard errors
#'@param R_rows Correlation matrix for rows of E
#'@param R_cols Correlation matrix for cols of E. Only one of R_row or R_cols can be given
#'@export
sim_bh3 <- function(true_L, true_F, true_Theta, S, R_rows, R_cols){
  M <- nrow(true_F)
  J <- nrow(true_L)
  K <- ncol(true_L)

  stopifnot(ncol(true_F) == K)
  stopifnot(nrow(true_Theta)==J & ncol(true_Theta) == M)
  stopifnot(nrow(S)==J & ncol(S) == M)
  if(!missing(R_rows) & !missing(R_cols)) stop("Please provide only one of R_rows or R_cols")
  if(missing(R_rows) & missing(R_cols)){
    rows <- TRUE
    R <- diag(rep(1, M))
  }else if(!missing(R_rows)){
    rows <- TRUE
    R <- R_rows
  }else{
    rows <- FALSE
    R <- R_cols
  }

  if(class(R)=="list" & rows==TRUE){
    blk_size <- purrr::map_dbl(R, function(b){nrow(b)})
    stopifnot(sum(blk_size)==J)
    nblock <- length(R)
    blk_start <- cumsum(c(1, blk_size[-nblock]))
    blk_stop <- blk_start + blk_size-1

    true_E <- purrr::map(seq_along(blk_size), function(i){
      apply(S[blk_start[i]:blk_stop[i],, drop=FALSE], 2, function(si){
        sigma <- diag(si) %*% R[[i]] %*% diag(si)
        MASS::mvrnorm(n=1, mu = rep(0, blk_size[i]), Sigma = sigma)
      })
    })
    true_E <- do.call(rbind, true_E)
  }
  else{
    if(!rows){
      stopifnot(nrow(R)==M & ncol(R) == M)
      true_E <- apply(S, 1, function(si){
        sigma <- diag(si) %*% R %*% diag(si)
        MASS::mvrnorm(n=1, mu = rep(0, M), Sigma = sigma)
      }) %>% t()
    }else{
      stopifnot(nrow(R)==J & ncol(R) ==J)
      true_E <- apply(S, 2, function(si){
        sigma <- diag(si) %*% R %*% diag(si)
        MASS::mvrnorm(n=1, mu = rep(0, J), Sigma = sigma)
      })
    }
  }
  beta_hat <- true_L %*% t(true_F) + true_Theta + true_E
  mats <- list(beta_hat = as.matrix(beta_hat), se_hat = S,
               true_L = true_L, true_F = true_F, true_E = true_E, true_Theta = true_Theta)
  return(mats)
}
