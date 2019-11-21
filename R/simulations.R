
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
