
#'@title Cheater ld-score matrix. Fast but higher variance.
#'@export
#'@keywords internal
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
#'@param Z_hat matrix of z-scores
#'@param ldscores vector of ldscores
#'@param ld_size Number of variants used to compute ldscores
#'@param N vector of sample sizes, length equal to number of columns of Z_hat
#'@param return_gencov If TRUE, the function will return genetic covariance in addition to residual covariance/correlation.
#'@param make_well_conditioned If TRUE, residual covariance will be projected to the nearest positive definite
#'matrix with condition number at least cond_num
#'@param cond_num Condition number used if make_well_conditioned = TRUE
#'@param blocks If not NULL, use jackknife to estimate the standard error of estimates.
#'@param ncores Number of cores to use for jackknifing
#'@return A list with some or all of the following elements. Se: Estimate of residual covariance Ve: Variance of Se.
#' Sg: genetic covariance Vg: Variance fo Sg. Rg: Genetic correlation, VRg: Variance of Rg.
#'@export
R_ldsc <- function(Z_hat,
                   ldscores,
                   ld_size,
                   N,
                   return_gencov = FALSE,
                   make_well_conditioned = TRUE,
                   cond_num = 1e5-1,
                   blocks = NULL,
                   ncores = 1){
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

  # h2 <- lapply(1:M, function(m){
  #   ix <- which(!is.na(Z_hat[,m]))
  #   snp_ldsc(ld_score = ldscores[ix], ld_size = ld_size, chi2 = Z_hat[ix,m]^2, sample_size = N[ix,m], blocks = NULL)
  # })
  res <- expand.grid(trait1 = 1:M, trait2 = 1:M) %>%
    filter(trait1 <= trait2)

  vals <- map2(res$trait1, res$trait2, function(i, j){
    #cat(i, j, "\n")
    # if(i == j){
    #   ix <- which(!is.na(Z_hat[,i]))
    #   h2 <- snp_ldsc(ld_score = ldscores[ix], ld_size = ld_size, chi2 = Z_hat[ix,i]^2, sample_size = N[ix,i], blocks = blocks)
    #   return(c(h2, 1))
    # }
    ix <- which(!is.na(Z_hat[,i]) & !is.na(Z_hat[,j]))
    rg <- ldsc_rg(ld_score = ldscores[ix], ld_size = ld_size,
            z1 = Z_hat[ix,i], z2 = Z_hat[ix,j],
            #h2_1 = h2[[i]], h2_2 = h2[[j]],
            h2_1 = NULL, h2_2 = NULL,
            sample_size_1 = N[ix,i], sample_size_2 = N[ix,j],
            blocks = blocks,
            ncores = ncores)
    rg
  })

  val_int <- map(vals, 1) %>% unlist()

  res$value <- val_int
  res_copy <- filter(res, trait1 != trait2) %>%
    rename(n1c = trait2, n2c = trait1) %>%
    rename(trait1 = n1c, trait2 = n2c)

  Se <- bind_rows(res, res_copy)  %>%
    reshape2::dcast(trait1 ~ trait2, value.var = "value")
  Se <- as.matrix(Se[,-1])
  colnames(Se) <- rownames(Se) <- NULL

  if(make_well_conditioned){
    #Se <- condition(Se, cond_num)
    Se <- Matrix::nearPD(Se, corr = FALSE, keepDiag = TRUE,
                         posd.tol = 1/cond_num)$mat
    Se <- as.matrix(Se)
  }

  if(!return_gencov & is.null(blocks)){
    return(list("Se" = Se))
  }

  if(!is.null(blocks)){
    val_resid_ve <- map(vals, 2) %>% unlist()
    val_gencov <- map(vals, 3) %>% unlist()
    val_gencor <- map(vals, 11) %>% unlist()
    val_gencov_ve <- map(vals, 4) %>% unlist()
    val_gencor_ve <- map(vals, 12) %>% unlist()
  }else{
    val_gencov <- map(vals, 2) %>% unlist()
    val_gencor <- map(vals, 7) %>% unlist()
  }

  ret <- list("Se" = Se)

  if(return_gencov){
    ## genetic covariance matrix
    res$value <- val_gencov
    res_copy <- filter(res, trait1 != trait2) %>%
      rename(n1c = trait2, n2c = trait1) %>%
      rename(trait1 = n1c, trait2 = n2c)

    Sg <- bind_rows(res, res_copy)  %>%
      reshape2::dcast(trait1 ~ trait2, value.var = "value")
    ret$Sg <- as.matrix(Sg[,-1])
    colnames(ret$Sg) <- rownames(ret$Sg) <- NULL


    ## genetic correlation matrix
    res$value <- val_gencor
    res_copy <- filter(res, trait1 != trait2) %>%
      rename(n1c = trait2, n2c = trait1) %>%
      rename(trait1 = n1c, trait2 = n2c)

    Rg <- bind_rows(res, res_copy)  %>%
      reshape2::dcast(trait1 ~ trait2, value.var = "value")
    ret$Rg <- as.matrix(Rg[,-1])
    colnames(ret$Rg) <- rownames(ret$Rg) <- NULL
  }
  if(!is.null(blocks)){
    res$value <- val_resid_ve^2
    res_copy <- filter(res, trait1 != trait2) %>%
      rename(n1c = trait2, n2c = trait1) %>%
      rename(trait1 = n1c, trait2 = n2c)

    Ve <- bind_rows(res, res_copy)  %>%
      reshape2::dcast(trait1 ~ trait2, value.var = "value")
    ret$Ve <- as.matrix(Ve[,-1])
    colnames(ret$Ve) <- rownames(ret$Ve) <- NULL
    if(return_gencov){
      ## genetic covariance matrix
      res$value <- val_gencov_ve^2
      res_copy <- filter(res, trait1 != trait2) %>%
        rename(n1c = trait2, n2c = trait1) %>%
        rename(trait1 = n1c, trait2 = n2c)

      Vg <- bind_rows(res, res_copy)  %>%
        reshape2::dcast(trait1 ~ trait2, value.var = "value")
      ret$Vg <- as.matrix(Vg[,-1])
      colnames(ret$Vg) <- rownames(ret$Vg) <- NULL


      ## genetic correlation matrix
      res$value <- val_gencor_ve^2
      res_copy <- filter(res, trait1 != trait2) %>%
        rename(n1c = trait2, n2c = trait1) %>%
        rename(trait1 = n1c, trait2 = n2c)

      VRg <- bind_rows(res, res_copy)  %>%
        reshape2::dcast(trait1 ~ trait2, value.var = "value")
      ret$VRg <- as.matrix(VRg[,-1])
      colnames(ret$VRg) <- rownames(ret$VRg) <- NULL
    }
  }
  return(ret)
}


#'@title Calculate matrix of error correlations using p-value threshold method.
#'@param B_hat Matrix of effect size estimates
#'@param S_hat Matrix of standard errors
#'@param p_val_thresh P-value threshold
#'@param return_cor Return matrix as correlation matrix (apply cov2cor), recommended
#'@param max_well_conditioned Project matrix to nearest well conditioned positive definite
#'matrix using Matrix::nearPD
#'@param cond_num Condition number if using make_well_conditioned
#'@return A traits by traits matrix of nuisance correlations.
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

  if(make_well_conditioned & !return_cor){
    #R <- condition(R, cond_num)
    R <- Matrix::nearPD(R, corr = FALSE, keepDiag = TRUE,
                         posd.tol = 1/cond_num)$mat
    R <- as.matrix(R)
  }else if(make_well_conditioned){
    R <- Matrix::nearPD(R, corr = TRUE,
                         posd.tol = 1/cond_num)$mat
    R <- as.matrix(R)
  }else if(return_cor){
    R <- cov2cor(R)
  }
  return(R)
}

