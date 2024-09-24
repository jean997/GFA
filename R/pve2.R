
# Should return a matrix that is traits by factors
pve2 <- function(gfa_fit){

  ef2 <- flashier:::get.EF2(gfa_fit$fit$flash_fit)
  k <- ncol(ef2[[1]])
  p <- nrow(ef2[[2]])
  n <- nrow(ef2[[1]])

  if(gfa_fit$method == "fixed_factors"){
    fixed_ix <- which(gfa_fit$fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    est_ix <- (1:gfa_fit$fit$n_factors)[-fixed_ix]
    #cat("est_ix:", est_ix, "\n")
    # relative variance explained per trait/factor
    sj <- sapply(est_ix, function(kk){
      Vk <- flashier:::lowrank.expand(list(ef2[[1]][,kk,drop = F], ef2[[2]][,kk, drop = F])) ## total effects on z-score scale
      ## standardized effect scale total variance
      colSums(Vk)/(gfa_fit$scale^2)
    }) |> matrix(nrow = p, byrow = FALSE)

    # variance from fixed factors (this is part of error)
    if(length(fixed_ix) == 0){
      fj_sums <- rep(0, p)
    }else{
      fj_sums <- sapply(fixed_ix, function(kk){
        Vk <- flashier:::lowrank.expand(list(ef2[[1]][,kk,drop = F], ef2[[2]][,kk, drop = F]))
        ## standardized effect scale total variance
        colSums(Vk)/(gfa_fit$scale^2)
      })|> matrix(nrow = p, byrow = FALSE) |> rowSums()
    }

    tau <- flashier:::get.tau(gfa_fit$fit$flash_fit)
    if(inherits(tau, "matrix")){
      stopifnot(nrow(tau) == n & ncol(tau) == p) # should never fail
      S2 <- gfa_fit$fit$flash_fit$given.S2
      est_tau <- 1/((1/tau) - S2)
    }else{
      stopifnot(length(tau) == p) # this should never fail
      # error component of tau
      fixed_tau <- flash_fit_get_fixed_tau(gfa_fit$fit$flash_fit)
      est_tau <- 1/((1/tau) - (1/fixed_tau))
    }
    error_var <- fj_sums + (n/fixed_tau)/(gfa_fit$scale^2)
    total_var <- rowSums(sj) + fj_sums + (n/tau)/(gfa_fit$scale^2)

    genet_var = rowSums(sj) + (n/est_tau)/(gfa_fit$scale^2)
    pve_j <- sj/genet_var
  }else if(gfa_fit$method == "random_effect"){
    est_ix <- 1:k
    # relative variance explained per trait/factor
    sj <- sapply(est_ix, function(kk){
      Vk <- flashier:::lowrank.expand(list(ef2[[1]][,kk,drop = F], ef2[[2]][,kk, drop = F])) ## total effects on z-score scale
      ## standardized effect scale total variance
      colSums(Vk)/(gfa_fit$scale^2)
    }) |> matrix(nrow = p, byrow = FALSE)
    est_tau <- flashier:::get.tau(gfa_fit$fit$flash_fit)
    error_var <- colSums(gfa_fit$fit$flash_fit$EB2)
    genet_var = rowSums(sj) + (n/est_tau)/(gfa_fit$scale^2)
    pve_j <- sj/genet_var
  }else if(gfa_fit$method == "noR"){
    est_ix <- 1:k
    # relative variance explained per trait/factor
    sj <- sapply(est_ix, function(kk){
      Vk <- flashier:::lowrank.expand(list(ef2[[1]][,kk,drop = F], ef2[[2]][,kk, drop = F])) ## total effects on z-score scale
      ## standardized effect scale total variance
      colSums(Vk)/(gfa_fit$scale^2)
    }) |> matrix(nrow = p, byrow = FALSE)

    tau <- flashier:::get.tau(gfa_fit$fit$flash_fit)
    if(inherits(tau, "matrix")){
      stopifnot(nrow(tau) == n & ncol(tau) == p) # should never fail
      S2 <- gfa_fit$fit$flash_fit$given.S2
      error_var <-  colSums(S2)/(gfa_fit$scale^2)
      tau_var <- colSums(1/tau)
      gvar_tau <- tau_var - error_var
    }else{
      stopifnot(length(tau) == p) # this should never fail
      # error component of tau
      fixed_tau <- flash_fit_get_fixed_tau(gfa_fit$fit$flash_fit)
      est_tau <- 1/((1/tau) - (1/fixed_tau))
      error_var <-  (n/fixed_tau)/(gfa_fit$scale^2)
      gvar_tau <- (n/est_tau)/(gfa_fit$scale^2)
    }
    #total_var <- rowSums(sj) + (n/tau)/(gfa_fit$scale^2)

    genet_var = rowSums(sj) + gvar_tau
    pve_j <- sj/genet_var
  }
  return(list(genet_var = genet_var, pve = pve_j))
}




