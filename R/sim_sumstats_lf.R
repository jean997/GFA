
#See simulating_data.Rmd for explanation
#'@title Simulate summary stats
#'@description Simulate summary statistcs for fully overlapping GWAS with no LD
#'@param F_mat factor matrix M by K
#'@param N GWAS sample size
#'@param J Total number of SNPs to generate
#'@param h_2_trait Heritability of each trait. Length M vector.
#'@param omega Proportion of SNP heritability mediated by factors. Length M vector.
#'@param h_2_factor Heritability of each factor. Length K vector.
#'@param pi_L Proportion of non-zero elements in L_k. Length K factor
#'@param pi_theta Proportion of non-zero elements in theta. Scalar.
#'@param R_E Environmental trait correlation not mediated by factors. M by M pd matrix
#'@param R_LD_list List of eigen decompositions of LD correlation matrices, may be missing.
#'@param relative_pve Relative pve contributed by each factor. Length K.
#'@param g_F Function from which non-zero elements of F are generated
#'@param pi_F Propotion of non-zero elements of F.
#'@details
#'
#'note: have removed relative_pve so below is not correct. will update,
#'omega, h2_trait, and relative_pve constrain the rows and columns of F. The sum of squared elements in row m must equal
#'omega[m]*h2_trait[m]. relative_pve is normalized to sum to sum(omega*h2_trait). After this normalization the sum of
#'squared elements in column k must equal relative_pve[k]. The function `scale_F` iteratively rescales rows and columns
#'until they have the desired sums. This allows us to preserve the input sparsity structure and approximate distribution
#'magnitudes and still achieve the desired row and column squared sums.
#'
#'If F_mat is provided it is rescaled as described. In this case relative_pve is optional.
#'
#'If F_mat is not provided, it will be generated using the `generate_F2` function.
#'In this case g_F and nz_factor must be provided. scale_factor is an optional argument that can
#'be use to adjust the relative scalings of the factors. All of the elements
#'in F are generated iid from a mixture of a point mass at 0 and g_F. The matrix is then
#'rescaled according to the constraints.
#'
#'
#'With this setup it is possible to specify a setting that is impossible. Usually this occurs if the heritability
#'of the factors is low but the heritability of the traits is high leading to a contradiction. Right now the function
#'will just return an error if that happens.
#'
#'
#'@export
sim_sumstats_lf <- function(F_mat, N, J, h_2_trait, omega, h_2_factor, pi_L, pi_theta,
                            R_E, R_LD,
                            g_F, nz_factor, scale_factor, add=FALSE,
                            overlap_prop =1){

  #Check parameters
  if(!missing(F_mat)){
    stopifnot("matrix" %in% class(F_mat))
    M <- nrow(F_mat)
    K <- ncol(F_mat)
  }else{
    if(missing(g_F) | missing(nz_factor) ){
      stop("If F_mat is missing please supply g_F and nz_factor")
    }
    K <- length(nz_factor)
    M <- length(h_2_trait)
    if(missing(scale_factor)) scale_factor <- rep(1, K)
    cat("Will generate L and F with ", J, " markers, ", M, " traits, and ", K, "factors.\n")
  }
  stopifnot(length(h_2_trait) == M)
  stopifnot(all(h_2_trait <= 1 & h_2_trait >= 0))
  stopifnot(length(h_2_factor) == K)
  stopifnot(all(h_2_factor <= 1 & h_2_factor >= 0))
  stopifnot(length(omega) == M)
  stopifnot(all(omega >= 0 & omega <= 1))
  stopifnot(length(pi_L) == K)
  stopifnot(all(pi_L <= 1 & pi_L > 0))
  stopifnot(length(pi_theta) == 1 )
  stopifnot(pi_theta >=0 & pi_theta <=1)
  stopifnot(nrow(R_E) == M & ncol(R_E) == M)
  stopifnot(Matrix::isSymmetric(R_E))
  stopifnot(length(scale_factor)==K)
  if(!missing(R_LD)) cat(R_LD, "\n")
  R_E_eig <- eigen(R_E)
  stopifnot(all(R_E_eig$values >= 0))
  if(any(omega < 1) & pi_theta == 0){stop("Theta contributes non-zero heritability so pi_theta must be greater than 0.")}

  if(!missing(F_mat)){
    if(any(rowSums(F_mat == 0) == K & omega >0)){
      stop("One row of F is all zero but corresponds to non-zero omega\n")
    }
  }
  #Re-scale F or generate it if it is missing
  if(missing(F_mat)){
    #relative_pve <- relative_pve*sum(omega*h_2_trait)/sum(relative_pve)
    #F_mat <- generate_F(percent_zero = 1-pi_F, square_row_sums = omega*h_2_trait,
    #                    square_col_sums = relative_pve, rfunc = g_F)
    F_mat <- generate_F2(non_zero_by_factor = nz_factor,
                         square_row_sums = omega*h_2_trait,
                         scale_by_factor = scale_factor,
                         rfunc = g_F, add=add)
    if(ncol(F_mat) > K){
      nextra <- ncol(F_mat)-K
      h_2_factor <- c(h_2_factor, rep(1, nextra))
      pi_L <- c(pi_L, rep(pi_theta, nextra))
    }
    if(any(rowSums(F_mat^2) == 0)){
      ix <- which(rowSums(F_mat^2)==0)
      omega[ix] <- 0
    }
  }else{
    #Re scale rows of F
    F_mat <- F_mat*sqrt(omega*h_2_trait/rowSums(F_mat^2))
  }


  #Generate theta
  sigma_theta <- sqrt( (1/(pi_theta*J))*(1-omega)*h_2_trait)
  sigma_theta[omega == 1] <- 0

  theta <- purrr::map(sigma_theta, function(s){
    t <- rbinom(n=J, size=1, prob = pi_theta)
    n <- sum(t==1)
    t[t==1] <- rnorm(n=n, mean=0, sd = s)
    return(t)
  }) %>% do.call(cbind, .)



  #Generate L
  sigma_L <- (1/(pi_L*J))
  L_mat <- purrr::map(pi_L, function(p){
    l <- rbinom(n=J, size=1, prob = p)
    n <- sum(l==1)
    l[l==1] <- rnorm(n=n, mean=0, sd = sqrt(1/(p*J)))
    return(l)
  }) %>% do.call(cbind, .)

  #Compute Beta
  beta = L_mat %*% t(F_mat) + theta

  #Compute row covariance
  Sigma_G <- F_mat %*% t(F_mat) + J*diag(pi_theta*sigma_theta^2)
  sigma_2_F <- (1-h_2_factor)/(h_2_factor)
  Sigma_FE <- F_mat %*% diag(sigma_2_F) %*% t(F_mat)
  sigma_E <- sqrt(1 - h_2_trait - diag(Sigma_FE))
  Sigma_E <- diag(sigma_E) %*% R_E %*% diag(sigma_E)
  Sigma <- (1/N)*(Sigma_G + Sigma_FE  + Sigma_E)
  Sigma_indep <- (1/N)*diag(M)
  Sigma <- overlap_prop*Sigma + (1-overlap_prop)*Sigma_indep

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

scale_F <- function(F_init, square_row_sums, square_col_sums, tol = 1e-5, max_rep = 100){
  stopifnot(all(square_row_sums >= 0))
  stopifnot(all(square_col_sums >= 0))
  stopifnot(abs(sum(square_row_sums) - sum(square_col_sums)) < tol)
  M <- length(square_row_sums)
  K <- length(square_col_sums)
  F_mat <- F_init
  r2 <- rowSums(F_mat^2)
  c2 <- colSums(F_mat^2)
  test <- max(abs(c(r2-square_row_sums, c2 - square_col_sums)))
  rep <- 1
  while(test > tol & rep <= max_rep){
    F_mat <- F_mat*sqrt(square_row_sums)/sqrt(r2)
    c2 <- colSums(F_mat^2)
    F_mat <- t(t(F_mat)*sqrt(square_col_sums)/sqrt(c2))
    r2 <- rowSums(F_mat^2)
    c2 <- colSums(F_mat^2)
    test <- max(abs(c(r2-square_row_sums, c2 - square_col_sums)))
    rep <- rep+1
    cat(rep, ": ", test," ", test > tol & rep <= max_rep, "\n")
  }
  return(F_mat)
}

generate_F <- function(percent_zero, square_row_sums, square_col_sums, rfunc = function(n){runif(n, -1, 1)},
                       tol = 1e-5, max_rep = 100){
  f <- function(n){
    x <- rbinom(n, 1, 1-percent_zero)
    x[x==1] <- rfunc(sum(x))
    return(x)
  }
  stopifnot(all(square_row_sums > 0))
  stopifnot(all(square_col_sums > 0))
  stopifnot(abs(sum(square_row_sums) - sum(square_col_sums)) < tol)
  M <- length(square_row_sums)
  K <- length(square_col_sums)

  init <- FALSE
  F_mat <- replicate(n=K, expr={f(M)})

  missing_ix <- which(rowSums(F_mat !=0)==0)
  square_row_sums <- square_row_sums[-missing_ix]

  F_mat <- scale_F(F_mat, square_row_sums, square_col_sums, tol, max_rep )

}

#'@export
generate_F2 <- function(non_zero_by_factor,
                        square_row_sums,
                        scale_by_factor,
                        rfunc = function(n){runif(n, -1, 1)},
                        add=FALSE){
  f <- function(n, nz){
    stopifnot(nz >= 1)
    x <- rep(0, n)
    ix <- sample(seq(n), size=nz, replace=F)
    x[ix] <- rfunc(nz)
    return(x)
  }
  stopifnot(all(square_row_sums > 0))

  M <- length(square_row_sums)
  K <- length(non_zero_by_factor)
  stopifnot(length(scale_by_factor)==K)

  F_mat <- sapply(seq(K), function(k){f(M, non_zero_by_factor[k])})
  F_mat <- t(t(F_mat)*scale_by_factor) # Multiply each column by scale factor
  if(any(rowSums(F_mat !=0)==0)){
    missing_ix <- which(rowSums(F_mat !=0)==0)
    if(add){
      # Add any missing traits
      F_add <- sapply(missing_ix, function(i){
                  x <- rep(0, M)
                  x[i] <- 1
                  return(x)
      })
      F_mat <- cbind(F_mat, F_add)
    }else{
      r2 <- rowSums(F_mat^2)
      F_mat <- F_mat*sqrt(square_row_sums)/sqrt(r2)
      F_mat[missing_ix,] <- 0
      return(F_mat)
    }
  }
  r2 <- rowSums(F_mat^2)
  F_mat <- F_mat*sqrt(square_row_sums)/sqrt(r2)
  return(F_mat)

}


