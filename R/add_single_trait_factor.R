
#'@export
add_single_trait_factor <- function(gfa_fit, ix){
  flash_fit <- gfa_fit[["fit"]][["flash_fit"]]

  est_ix <- 1:ncol(gfa_fit$F_hat)
  n <- length(est_ix)
  error_ix <- gfa_fit$error_ix

  if(any(est_ix %in% error_ix)){
    stop("Something is wrong with factor indexing.\n")
  }
  if(gfa_fit$num_single > 0){
    single_ix <- (n + 1):(n + gfa_fit$num_single)
    n <- n + gfa_fit$num_single
    single_traits <- which(rowSums(gfa_fit$F_hat_single) != 0)
    if(all(ix %in% single_traits)){
      warning("Requested single-trait factors are already present")
      return(gfa_fit)
    }
    ix <- ix[!ix %in% single_traits]
    if(any(single_ix %in% error_ix)){
      stop("Something is wrong with factor indexing.\n")
    }
  }
  N <- gfa_fit[["fit"]][["n_factors"]]
  if(!all(ix %in% 1:ncol(flash_fit[["Y"]]))){
    stop("ix should be between 1 and ", ncol(flash_fit[["Y"]]), "\n")
  }

  lft <- flashier:::lowrank.expand(list(flash_fit[["EF"]][[1]], flash_fit[["EF"]][[2]]))
  resid <- flash_fit[["Y"]] - lft
  stF <- matrix(0, nrow = ncol(resid), ncol = length(ix))
  for(i in seq_along(ix)){
    stF[ix[i], i] <- 1
  }

  fitn <- flash_fit %>%
    flash_factors_init(init = list(resid[,ix, drop = FALSE], stF),
                       ebnm_fn = list(gfa_fit$params$ebnm_fn_L, gfa_fit$params$ebnm_fn_F)) %>%
    #flash_factors_reorder(new_order) %>%
    flash_factors_fix(., kset = (N + 1):(N + length(ix)), which_dim = "factors") %>%
    flash_backfit()

  new_single_ix <- (N + 1):(N + length(ix))
  new_order <- c(est_ix, single_ix, new_single_ix, error_ix )

  fitn <- flash_factors_reorder(fitn, new_order) %>%
          flash_nullcheck
  fix.dim <- flashier:::get.fix.dim(fitn$flash_fit)
  if(length(fix.dim) == 0){
    fixed_ix <- rep(FALSE, nfactor)
    num_single_fixed <- 0
  }else{
    fixed_ix <- fix.dim %>% sapply(., function(x){
      if(is.null(x)) return(FALSE)
      if(x == 2) return(TRUE)
      return(FALSE)})
    fixed_F <- fitn$F_pm[,fixed_ix == TRUE]
    nf <- apply(fixed_F, 2, function(x){sum(x !=0)})
    num_single_fixed <- sum(nf == 1)
  }
  gfa_fit_new <- gfa_wrapup(fitn, method = gfa_fit[["method"]], scale = gfa_fit[["scale"]],
                            num_single_fixed = num_single_fixed, nullcheck = FALSE)
}
