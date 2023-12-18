
#'@title Cheater ld-score matrix. Fast but higher variance.
#'@export
R_ldsc_quick <- function(Z_hat, ldscores, weights = 1/ldscores,
                         make_well_conditioned = TRUE, cond_num = 1e5,
                         return_cor = TRUE){
  M <- ncol(Z_hat)
  stopifnot(class(ldscores) == "numeric")
  stopifnot(length(ldscores) == nrow(Z_hat))
  res <- expand.grid(trait1 = 1:M, trait2 = 1:M) %>%
         filter(trait1 <= trait2)
  vals <- map2_dbl(res$trait1, res$trait2, function(i, j){
          zz <- Z_hat[,i]*Z_hat[,j]
          f <- lm(zz ~ ldscores, weights = weights)
          return(f$coefficients[1])})
  res$value <- vals
  res_copy <- filter(res, trait1 != trait2) %>%
    rename(n1c = trait2, n2c = trait1) %>%
    rename(trait1 = n1c, trait2 = n2c)

  cov_mat <- bind_rows(res, res_copy)  %>%
    reshape2::dcast(trait1 ~ trait2, value.var = "value")
  Re <- as.matrix(cov_mat[,-1])


  if(make_well_conditioned){
    Re <- condition(Re, cond_num)
  }

  if(return_cor){
    Re <- cov2cor(Re)
  }

  return(Re)
}

#'@title Calculate matrix of error correlations using LD-score regression.
#'@export
R_ldsc <- function(Z_hat, ldscores, ld_size, N, return_gencov = FALSE,
                   make_well_conditioned = TRUE, cond_num = 1e5-1,
                   return_cor = TRUE){
  M <- ncol(Z_hat)
  J <- nrow(Z_hat)
  stopifnot(class(ldscores) == "numeric")
  stopifnot(length(ldscores) == J)
  if("matrix" %in% class(N)){
    stopifnot(nrow(N) == J & ncol(N) == M)
  }else if(class(N) == "numeric"){
    stopifnot(length(N) == M)
    N <- matrix(rep(N, each = J), nrow = J)
  }

  h2 <- lapply(1:M, function(m){
    snp_ldsc(ld_score = ldscores, ld_size = ld_size, chi2 = Z_hat[,m]^2, sample_size = N[,m], blocks = NULL)
  })
  res <- expand.grid(trait1 = 1:M, trait2 = 1:M) %>%
    filter(trait1 <= trait2)

  vals <- map2(res$trait1, res$trait2, function(i, j){
    if(i == j) return(h2[[i]])
    rg <- ldsc_rg(ld_score = ldscores, ld_size = ld_size,
            z1 = Z_hat[,i], z2 = Z_hat[,j],
            h2_1 = h2[[i]], h2_2 = h2[[j]],
            sample_size_1 = N[,i], sample_size_2 = N[,j],
            blocks = NULL)
    c(rg[["int"]], rg[["gencov"]], rg[["gencorr"]])
  })

  val_int <- map(vals, 1) %>% unlist()

  res$value <- val_int
  res_copy <- filter(res, trait1 != trait2) %>%
    rename(n1c = trait2, n2c = trait1) %>%
    rename(trait1 = n1c, trait2 = n2c)

  cov_mat <- bind_rows(res, res_copy)  %>%
    reshape2::dcast(trait1 ~ trait2, value.var = "value")
  Re <- as.matrix(cov_mat[,-1])
  colnames(Re) <- rownames(Re) <- NULL

  if(make_well_conditioned){
    Re <- condition(Re, cond_num)
  }

  if(return_cor){
    Re <- cov2cor(Re)
  }

  if(!return_gencov){return(Re)}

  res_gencov <- expand.grid(trait1 = 1:M, trait2 = 1:M) %>%
    filter(trait1 <= trait2)
  res_gencov$value <- map(vals, 2) %>% unlist()
  res_g_copy <- filter(res_gencov, trait1 != trait2) %>%
    rename(n1c = trait2, n2c = trait1) %>%
    rename(trait1 = n1c, trait2 = n2c)

  gcov_mat <- bind_rows(res_gencov, res_g_copy)  %>%
    reshape2::dcast(trait1 ~ trait2, value.var = "value")
  Rg <- as.matrix(gcov_mat[,-1])

  res_gencor <- expand.grid(trait1 = 1:M, trait2 = 1:M) %>%
    filter(trait1 <= trait2)
  res_gencor$value <- map(vals, 3) %>% unlist()
  res_g_copy <- filter(res_gencor, trait1 != trait2) %>%
    rename(n1c = trait2, n2c = trait1) %>%
    rename(trait1 = n1c, trait2 = n2c)

  gcor_mat <- bind_rows(res_gencor, res_g_copy)  %>%
    reshape2::dcast(trait1 ~ trait2, value.var = "value")
  Rgcor <- as.matrix(gcov_mat[,-1])

  return(list("Re" = Re, "Rg" = Rg, "Rgcor" = Rgcor))

}


#'@title Calculate matrix of error correlations using p-value threshold method.
#'@export
R_pt <- function(B_hat, S_hat, p_val_thresh = 0.05, return_cor = TRUE,
                 make_well_conditioned = TRUE, cond_num = 1e5){
  ntrait <- ncol(B_hat)
  nvar <- nrow(B_hat)
  stopifnot(nrow(S_hat) == nvar & ncol(S_hat)==ntrait)

  keep <- 2*pnorm(-abs(B_hat/S_hat)) > p_val_thresh

  A <- B_hat/S_hat

  R_df <- expand.grid(T1 = seq(ntrait), T2 = seq(ntrait)) %>%
    filter(T1 <= T2)
  r <- map2_dbl(R_df$T1, R_df$T2, function(i, j){
    if(i == j) return(1)
    ix <- which(keep[,i]==TRUE & keep[,j]==TRUE)
    cor(A[ix,i], A[ix,j])
  } )
  R_df$r <-r
  R <- reshape2::acast(R_df, T1~T2, value.var = "r")
  R_lower <- reshape2::acast(R_df, T2~T1, value.var = "r")
  R[lower.tri(R)] <- R_lower[lower.tri(R)]

  if(make_well_conditioned){
    R <- condition(R, cond_num)
  }

  if(return_cor){
    R <- cov2cor(R)
  }

  return(R)
}

#'@title Project matrix to nearest well conditioned positive definite matrix.
#'@param R Matrix
#'@param cond_num Maximum allowable condition number (max(eigenvalue)/min(eigenvalue))
#'@export
condition <- function(R, cond_num = 1e5, corr = FALSE){
  eig_R <- eigen(R)
  vals <- eig_R$values
  illcond <- any(vals < 0) | max(vals)/min(vals) > cond_num
  if(illcond){
    warning(paste0(deparse(substitute(R)), " is either not positive definite or is illconditioned. Projecting to nearest well conditioned matrix."))
    x <- (cond_num*min(vals) - max(vals))/(1-cond_num)
    eig_R$values <- eig_R$values + x
    R <- with(eig_R, tcrossprod(vectors, tcrossprod(vectors, diag(values))))
  }
  if(corr){
    R <- cov2cor(R)
  }
  return(R)
}

