gfa_singletrait_check <- function(fit, check_thresh = 0.9, params){

  nfactor <- ncol(fit$F_pm)

  fix_ix <- my_flash_get_fixed_idx(fit)
  num_error_fixed <- length(fix_ix$factors)
  error_ix <- fix_ix$factors
  num_single_fixed <- 0
  num_est <- nfactor - num_error_fixed

  checked_trait <- c()

  done <- FALSE # stop if there are no single trait candidates or none were replaced
  while(!done){
    D <- fit$F_pm
    nfactor <- ncol(D)
    L <- fit$L_pm


    fix_ix <- my_flash_get_fixed_idx(fit)
    #fix.dim <- flashier:::get.fix.dim(fit$flash_fit)
    is_fixed <- rep(FALSE, nfactor)
    if(length(fix_ix$factors) > 0){
      is_fixed[fix_ix$factors] <- TRUE
    }

    Dn <- norm_cols(D)$A
    col_max <- apply(abs(Dn), 2, max)

    single_trait_index <- which(col_max > check_thresh & !is_fixed)
    single_traits <- sapply(single_trait_index, function(i){
      which.max(abs(Dn[,i]))
    })
    if(length(single_trait_index) == 0){
      done <- TRUE
    }else{
      new_single_fixed_ix <- c()
      for(i in single_trait_index){
        cat("Checking factor ", i, "\n")
        myfactor <- D[,i, drop = FALSE]
        myloadings <- L[, i, drop = FALSE]
        altfactor <- myfactor
        k <- which.max(abs(myfactor))
        altfactor[-k] <- 0
        fitn <- flash_factors_remove(fit, i) %>%
            flash_factors_init(init = list(myloadings, altfactor),
                               ebnm_fn = list(params$ebnm_fn_L, params$ebnm_fn_F)) %>%
            flash_factors_fix(., kset = nfactor, which_dim = "factors") %>%
            flash_backfit()

        ## arrange back into original order
        if(i < nfactor){
          new_order <- c(1:(i-1), nfactor, i:(nfactor-1))
          fitn <- flash_factors_reorder(fitn, new_order)
        }
        if(fitn$elbo > fit$elbo){
          message(paste0("Replacing factor ", i , " with single trait factor"))
          fit <- fitn
          new_single_fixed_ix <- c(new_single_fixed_ix, i)
        }
      }
      ## reorder after checking all so that single trait factors are at the end
      if(length(new_single_fixed_ix) > 0){
        new_order <- (1:num_est)
        new_order <- new_order[!new_order %in% new_single_fixed_ix]
        new_order <- c(new_order, new_single_fixed_ix, error_ix)
        fit <- flash_factors_reorder(fit, new_order)
        num_single_fixed <- num_single_fixed + length(new_single_fixed_ix)
      }else{
        done <- TRUE
      }
    }
  }
  return(list(fit = fit,  n_est = num_est - num_single_fixed, n_single = num_single_fixed, n_error = num_error_fixed))
}


