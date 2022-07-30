#'@export
R_ldsc <- function(Z_hat, ldscores, weights){
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
  return(Re)
}

#'@export
R_ldsc2 <- function(Z_hat, ldscores, ld_size, N){
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
  res <- expand.grid(trait1 = 1:M, trait2 = 1:M) %>%
    filter(trait1 <= trait2)
  vals <- map2_dbl(res$trait1, res$trait2, function(i, j){
    rg <- bigsnpr::snp_ldsc_rg(ld_score = ldscores, ld_size = ld_size,
                               z1 = Z_hat[,i], z2 = Z_hat[,j],
                               sample_size_1 = N[,i], sample_size_2 = N[,j],
                               blocks = NULL)
    return(rg[["int"]])
  })

  res$value <- vals
  res_copy <- filter(res, trait1 != trait2) %>%
    rename(n1c = trait2, n2c = trait1) %>%
    rename(trait1 = n1c, trait2 = n2c)

  cov_mat <- bind_rows(res, res_copy)  %>%
    reshape2::dcast(trait1 ~ trait2, value.var = "value")
  Re <- as.matrix(cov_mat[,-1])
  return(Re)
}
