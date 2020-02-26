
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
