
#See simulating_data.Rmd for explanation
#'@title Simulate summary stats
#'@description Simulate summary statistcs for fully overlapping GWAS with no LD
#'@param F_mat factor matrix M by K
#'@param N GWAS sample size
#'@param J Total number of SNPs to generate
#'@param h_2_trait Heritability of each trait. Length M vector.
#'@param omega Proportion of SNP heritability mediated by factors. Length M vector.
#'@param h_2_factor Heritability of each factor. Length K vector.
#'@param pi_L SNP sparsity of each factor. Length K factor
#'@param pi_theta SNP sparsity of theta. Scalar.
#'@param R_E Environmental trait correlation not mediated by factors. M by M pd matrix
#'@param R_LD_list List of eigen decompositions of LD correlation matrices, may be missing.
#'@export
sim_sumstats_lf <- function(F_mat, N, J, h_2_trait, omega, h_2_factor, pi_L, pi_theta,
                            R_E, R_LD){

  #Check parameters
  stopifnot("matrix" %in% class(F_mat))
  M <- nrow(F_mat)
  K <- ncol(F_mat)
  stopifnot(length(h_2_trait) == M)
  stopifnot(all(h_2_trait < 1 & h_2_trait >= 0))
  stopifnot(length(h_2_factor) == K)
  stopifnot(all(h_2_factor <= 1 & h_2_factor >= 0))
  stopifnot(length(omega) == M)
  stopifnot(all(omega >= 0 & omega <= 1))
  stopifnot(length(pi_L) == K)
  stopifnot(all(pi_L < 1 & pi_L > 0))
  stopifnot(length(pi_theta) == 1 )
  stopifnot(pi_theta >=0 & pi_theta <=1)
  stopifnot(nrow(R_E) == M & ncol(R_E) == M)
  stopifnot(Matrix::isSymmetric(R_E))
  R_E_eig <- eigen(R_E)
  stopifnot(all(R_E_eig$values >= 0))

  #if(!missing(R_LD)){
  #  l <- sapply(R_LD, function(e){length(e$values)})
  #  if(sum(l) != J){
  #    stop("LD information supplied for ", sum(l), " snps, not ", J, ".")
  #  }
  #}

  if(any(rowSums(F_mat == 0) == K & omega >0)){
    stop("One row of F is all zero but corresponds to non-zero omega\n")
  }


  #Re-scale F
  F_rowsums <- omega*h_2_trait
  x <- rowSums(F_mat^2)
  F_mat <- F_mat*sqrt(F_rowsums/x)

  #Generate L
  sigma_L <- (1/(pi_L*J))
  L_mat <- purrr::map_dfc(pi_L, function(p){
    l <- rbinom(n=J, size=1, prob = p)
    n <- sum(l==1)
    l[l==1] <- rnorm(n=n, mean=0, sd = sqrt(1/(p*J)))
    return(l)
  }) %>% as.matrix()

  #Generate theta
  sigma_theta <- sqrt( (1/(pi_theta*J))*(1-omega)*h_2_trait)
  theta <- purrr::map_dfc(sigma_theta, function(s){
    t <- rbinom(n=J, size=1, prob = pi_theta)
    n <- sum(t==1)
    t[t==1] <- rnorm(n=n, mean=0, sd = s)
    return(t)
  }) %>% as.matrix()

  #Compute Beta
  beta = L_mat %*% t(F_mat) + theta

  #Compute row covariance
  Sigma_G <- F_mat %*% t(F_mat) + J*diag(pi_theta*sigma_theta^2)
  sigma_2_F <- (1-h_2_factor)/(h_2_factor)
  Sigma_FE <- F_mat %*% diag(sigma_2_F^2) %*% t(F_mat)
  sigma_E <- sqrt(1 - h_2_trait - diag(Sigma_FE))
  Sigma_E <- diag(sigma_E) %*% R_E %*% diag(sigma_E)
  Sigma <- (1/N)*(Sigma_G + Sigma_FE  + Sigma_E)

  #Generate sampling error
  E <- MASS::mvrnorm(n=J, mu = rep(0, M), Sigma = Sigma)

  #Generate summary statistics
  if(missing(R_LD)){
    beta_hat <- beta + E
    se_beta_hat <- matrix(rep(sqrt(diag(Sigma)), J), byrow=T, nrow=J)

    ret <- list(beta_hat =beta_hat, se_beta_hat = se_beta_hat,
                L_mat = L_mat, F_mat = F_mat, theta = theta,
                R_E = R_E, Sigma=Sigma)
    return(ret)
  }

  #If LD, introduce correlation

  #Convert to z-scores
  S <- matrix(rep(sqrt(diag(Sigma)), J), byrow=T, nrow=J)
  E_Z <- E/S
  Z <- beta/S

  #Figure out how much/how many replicates of supplied LD we need
  # This is the fidliest bit
  l <- sapply(R_LD, function(e){length(e$values)})
  nblock <- length(R_LD)
  full_reps <- floor(J/sum(l))
  remainder <- J  - full_reps*sum(l)
  blocks_rem <- max(which(cumsum(l) <= remainder))
  final_remainder <- remainder-cumsum(l)[blocks_rem]

  last_block <- with(R_LD[[blocks_rem + 1]], (vectors %*% diag(values) %*% vectors)[1:final_remainder, 1:final_remainder])
  R_LD[[nblock + 1]] <- eigen(last_block)
  block_index <- c(rep(seq(nblock), full_reps), 1:blocks_rem, nblock + 1)
  l <- c(l, final_remainder)[block_index]
  start_ix <- cumsum(c(1, l[-length(l)]))
  end_ix <- start_ix + l-1

  E_LD_Z <- lapply(seq_along(block_index), function(i){
    with(R_LD[[block_index[i]]], vectors %*% sqrt(diag(values)) %*% E_Z[start_ix[i]:end_ix[i], ])
  }) %>% do.call( rbind, .)

  Z_hat <- Z + E_LD_Z

  beta_hat <- Z_hat*S
  E_LD <- E_LD_Z*S

  ret <- list(beta_hat =beta_hat, se_beta_hat = S,
              L_mat = L_mat, F_mat = F_mat, theta = theta,
              R_E = R_E, Sigma=Sigma)
  return(ret)
}
