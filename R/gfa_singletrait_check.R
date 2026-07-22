gfa_singletrait_check <- function(fit, check_thresh = 0.9, params){

  nfactor <- ncol(fit$F_pm)

  fix.dim <- flashier:::get.fix.dim(fit$flash_fit)
  if(length(fix.dim) == 0){
    fixed_ix <- rep(FALSE, nfactor)
  }else{
    fixed_ix <- fix.dim %>% sapply(., function(x){
      if(is.null(x)) return(FALSE)
      if(x == 2) return(TRUE)
      return(FALSE)})
  }
  num_error_fixed <- sum(fixed_ix)
  num_single_fixed <- 0
  num_est <- nfactor - num_error_fixed

  checked_trait <- c()
  upd_this_round <- FALSE
  done <- FALSE
  while(!done){
    D <- fit$F_pm
    nfactor <- ncol(D)
    L <- fit$L_pm

    fix.dim <- flashier:::get.fix.dim(fit$flash_fit)
    if(length(fix.dim) == 0){
      fixed_ix <- rep(FALSE, nfactor)
    }else{
      fixed_ix <- fix.dim %>% sapply(., function(x){
         if(is.null(x)) return(FALSE)
         if(x == 2) return(TRUE)
         return(FALSE)})
    }

    Dn <- norm_cols(D)$A
    col_max <- apply(abs(Dn), 2, max)

    single_trait_index <- which(col_max > check_thresh & !fixed_ix)
    single_traits <- sapply(single_trait_index, function(i){
      which.max(abs(Dn[,i]))
    })
    single_trait_index <- single_trait_index[!single_traits %in% checked_trait]
    if(length(single_trait_index) == 0){
      if(upd_this_round){
        checked_trait <- c()
        upd_this_round <- FALSE
      }else{
        done <- TRUE
      }
    }else{
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
            #flash_factors_reorder(new_order) %>%
            flash_factors_fix(., kset = nfactor, which_dim = "factors") %>%
            flash_backfit()



        #class(fitn$flash_fit$EF) <-  c("lowrank", "list")
        #class(fitn$flash_fit$EF2) <-  c("lowrank", "list")

        if(fitn$elbo > fit$elbo){
          message(paste0("Replacing factor ", i , " with single trait factor"))
          num_est <- num_est -1
          num_single_fixed <- num_single_fixed + 1

          new_order <- seq(nfactor)
          if(num_error_fixed > 0){
            new_order[(num_est+ 2):nfactor] <- (num_est + 1):(nfactor-1)
            new_order[num_est + 1] <- nfactor
          }
          fitn <- flash_factors_reorder(fitn, new_order)

          fit <- fitn
          upd_this_round <- TRUE
        }
        checked_trait <- c(checked_trait, k)

        fix.dim <- flashier:::get.fix.dim(fit$flash_fit)
        if(length(fix.dim) == 0){
          fixed_ix <- rep(FALSE, nfactor)
        }else{
          fixed_ix <- fix.dim %>% sapply(., function(x){
            if(is.null(x)) return(FALSE)
            if(x == 2) return(TRUE)
            return(FALSE)})
        }


      }
    }
  }
  fit$num_single_fixed <- num_single_fixed
  fit$num_error_fixed <- num_error_fixed
  return(fit)
}


