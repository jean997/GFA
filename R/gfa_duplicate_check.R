#'@export
gfa_duplicate_check <- function(fit, dim = 2, check_thresh = 0.5){
  if(!dim %in% c(1, 2)) stop("dim must be 1 or 2 in gfa_duplicate_check.\n")
  done <- FALSE
  while(!done){
    cat("begin while\n")
    if(dim == 2){
      D = fit$F_pm
    }else if(dim == 1){
      D = fit$L_pm
    }
    if(length(fit$flash_fit$fix.dim) == 0){
      fixed_ix <- rep(FALSE, ncol(D))
    }else{
      fixed_ix <- fit$flash_fit$fix.dim %>% sapply(., function(x){
          if(is.null(x)) return(FALSE)
          if(x == dim) return(TRUE)
          return(FALSE)})
    }

    if(any(fixed_ix)){
      fixed_ix <- which(fixed_ix)
      original_ix <- seq(ncol(D))
      D <- D[,-fixed_ix, drop=FALSE]
      non_fixed_ix = original_ix[-fixed_ix]
    }else{
      non_fixed_ix = seq(ncol(D))
    }

    Dn <- norm_cols(D)

    d <- t(Dn$A) %*% Dn$A
    diag(d) <- 0
    if(any(abs(d) > check_thresh)){
      dups <- reshape2::melt(abs(d)) %>%
              filter(value > check_thresh)

      dups$w1 <- Dn$w[dups$Var1]
      dups$w2 <- Dn$w[dups$Var2]
      dups <- filter(dups, w1 < w2) %>%
              arrange(-value)
      cat(dim(dups), "\n")
      print(dups)
      reset <- FALSE
      for(i in dups$Var1){
        cat(i, "\n")
        fit_rem <- flash_factors_remove(fit, kset = non_fixed_ix[i]) %>% flash_backfit()
        if(fit_rem$elbo > fit$elbo){
          cat("Removing factor ", non_fixed_ix[i], "\n")
          fit <- fit_rem
          if(any(fit$flash_fit$is.zero)){
            fit <- flash_nullcheck(fit, tol = -Inf)
          }
          reset <- TRUE
        }else{
          cat("Not removing factor ", non_fixed_ix[i], "\n")
        }
        if(reset) break
      }
      if(i == last(dups$Var1) & !reset){
        done <- TRUE
      }
    }else{
      done <- TRUE
    }
  }
  return(fit)
}
