
# Should return a matrix that is traits by factors
pve2 <- function(gfa_fit){

  ef2 <- flashier:::get.EF2(gfa_fit$fit$flash_fit)
  k <- ncol(ef2[[1]])
  p <- nrow(ef2[[2]])
  n <- nrow(ef2[[1]])

  if(gfa_fit$method == "fixed_factors"){
    fixed_ix <- which(gfa_fit$fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    est_ix <- which(!gfa_fit$fit$flash_fit$fix.dim %>% sapply(., function(x){!is.null(x)}))
    #cat("est_ix:", est_ix, "\n")
    # relative variance explained per trait/factor
    sj <- sapply(est_ix, function(kk){
      Vk <- flashier:::lowrank.expand(list(ef2[[1]][,kk,drop = F], ef2[[2]][,kk, drop = F])) ## total effects on z-score scale
      ## standardized effect scale total variance
      colSums(Vk)/(gfa_fit$scale^2)
    }) |> matrix(nrow = p, byrow = FALSE)

    # variance from fixed factors (this is part of error)
    fj <- sapply(fixed_ix, function(kk){
      Vk <- flashier:::lowrank.expand(list(ef2[[1]][,kk,drop = F], ef2[[2]][,kk, drop = F]))
      ## standardized effect scale total variance
      colSums(Vk)/(gfa_fit$scale^2)
    })|> matrix(nrow = p, byrow = FALSE)

    tau <- flashier:::get.tau(gfa_fit$fit$flash_fit)
    stopifnot(length(tau) == p) # this should never fail
    # error component of tau
    fixed_tau <- flash_fit_get_fixed_tau(gfa_fit$fit$flash_fit)
    est_tau <- 1/((1/tau) - (1/fixed_tau))

    error_var <- rowSums(fj) + (n/fixed_tau)/(gfa_fit$scale^2)
    total_var <- rowSums(sj) + rowSums(fj) + (n/tau)/(gfa_fit$scale^2) ## this should be close to 1?

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
    error_var <- colSums(gfa_fit$fit$flash_fit$EB2)
    genet_var <- rowSums(sj) + (n/flashier:::get.tau(gfa_fit$fit$flash_fit))
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
    stopifnot(length(tau) == p) # this should never fail
    # error component of tau
    fixed_tau <- flash_fit_get_fixed_tau(gfa_fit$fit$flash_fit)
    est_tau <- 1/((1/tau) - (1/fixed_tau))

    error_var <-  (n/fixed_tau)/(gfa_fit$scale^2)
    total_var <- rowSums(sj) + (n/tau)/(gfa_fit$scale^2) ## this should be close to 1?

    genet_var = rowSums(sj) + (n/est_tau)/(gfa_fit$scale^2)
    pve_j <- sj/genet_var
  }
  return(list(geent_var = genet_var, error_var = error_var, pve = pve_j))
}




