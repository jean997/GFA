
# Should return a matrix that is traits by factors
pve2 <- function(gfa_fit){

  ef2 <- flashier:::get.EF2(gfa_fit$fit$flash_fit)
  k <- ncol(ef2[[1]])
  p <- nrow(ef2[[2]])
  n <- nrow(ef2[[1]])

  tau <- flashier:::get.tau(gfa_fit$fit$flash_fit)
  if(inherits(tau, "matrix")){
    stopifnot(nrow(tau) == n & ncol(tau) == p) # should never fail
    S2 <- gfa_fit$fit$flash_fit$given.S2
    gvar_theta <- (colSums(1/tau) - colSums(S2))/gfa_fit$scale^2
  }else{
    stopifnot(length(tau) == p) # this should never fail
    # error component of tau
    fixed_tau <- flash_fit_get_fixed_tau(gfa_fit$fit$flash_fit)
    est_tau <- 1/((1/tau) - (1/fixed_tau))
    #error_var <-  (n/fixed_tau)/(gfa_fit$scale^2)
    gvar_theta <- (n/est_tau)/(gfa_fit$scale^2)
  }

  est_ix <- 1:k
  if(gfa_fit$method == "fixed_factors"){
    fixed_ix <- which(gfa_fit$fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    if(length(fixed_ix) > 0){
      est_ix <- est_ix[-fixed_ix]
    }
  }

  sj <- sapply(est_ix, function(kk){
      Vk <- flashier:::lowrank.expand(list(ef2[[1]][,kk,drop = F], ef2[[2]][,kk, drop = F])) ## total effects on z-score scale
      ## standardized effect scale total variance
      colSums(Vk)/(gfa_fit$scale^2)
  }) |> matrix(nrow = p, byrow = FALSE)

  genet_var = rowSums(sj) + gvar_theta
  pve_j <- sj/genet_var

  return(list(genet_var = genet_var, pve = pve_j))
}

### PVE calculation notes:
#
## We fit
## Zhat = Bhatstd*S = L ( S F)^T + Theta S + E_s
## where S is diagonal matrix of scales.
## For continuous trait, the scale is sqrt(N)
## Bhatstd gives effect estimates on the sd(Y)/sd(G) scale.
##
## For the noR case, elements of E_s are independent standard normal.
## for the fixed factors case, E_s = Eprime + AW where W are fixed factors
##
## The genetic variance of trait j explained by factor k is
##  sum_i (L_{i,k} F_{j,k})^2. The posterior expected value of this is given by
# Vk <- flashier:::lowrank.expand(list(ef2[[1]][,k,drop = F], ef2[[2]][,k, drop = F])) ## total effects on z-score scale
# colSums(Vk)/(gfa_fit$scale^2) ## This has length p
##
## The genetic variance explained by theta is stored in the tau object.
## (1/tau) = colSums(theta_s^2)/n + colSums(Eprime^2)/n
## colSums(Eprime^2)/n = 1/fixed_tau (= 1 noR or lambda_min for fixed_factors)
## Therefore, the genetic variance explained by tau is
## n((1/tau) - (1/fixed_tau)) = colSums(theta_s^2)
## n((1/tau) - (1/fixed_tau))/scale^2 = colSums(theta^2)
## This is given in gvar_theta
## so total genetic variance for trait j is rowSums(sj) + gvar_theta

## The input sample has already been scaled to match L_hat and F_hat
pve_from_sample <- function(gfa_fit, samp){

  k <- ncol(samp[["L"]])
  p <- nrow(samp[["F"]])
  n <- nrow(samp[["L"]])

  tau <- flashier:::get.tau(gfa_fit$fit$flash_fit)
  if(inherits(tau, "matrix")){
    stopifnot(nrow(tau) == n & ncol(tau) == p) # should never fail
    S2 <- gfa_fit$fit$flash_fit$given.S2
    gvar_theta <- (colSums(1/tau) - colSums(S2))/gfa_fit$scale^2
  }else{
    stopifnot(length(tau) == p) # this should never fail
    # error component of tau
    fixed_tau <- flash_fit_get_fixed_tau(gfa_fit$fit$flash_fit)
    est_tau <- 1/((1/tau) - (1/fixed_tau))
    #error_var <-  (n/fixed_tau)/(gfa_fit$scale^2)
    gvar_theta <- (n/est_tau)/(gfa_fit$scale^2)
  }

  est_ix <- 1:k
  if(gfa_fit$method == "fixed_factors"){
    fixed_ix <- which(gfa_fit$fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    if(length(fixed_ix) > 0){
      est_ix <- est_ix[-fixed_ix]
    }
  }

  sj <- sapply(est_ix, function(kk){
    Vk <- flashier:::lowrank.expand(list((samp[["L"]][,kk,drop = F])^2, (samp[["F"]][,kk, drop = F])^2)) ## total effects on z-score scale
    ## standardized effect scale total variance
    colSums(Vk) #/(gfa_fit$scale^2) not needed because sample is already scaled
  }) |> matrix(nrow = p, byrow = FALSE)

  genet_var = rowSums(sj) + gvar_theta
  pve_j <- sj/genet_var


  return(list(genet_var = genet_var, pve = pve_j))
}
