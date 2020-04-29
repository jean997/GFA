

#'@title Simulate Individual Level Data
#'@param nind Number of individuals
#'@param nvar Number of variants
#'@param F_mat Trait loadings (traits by factor)
#'@param pi_lodaings Proportion of non-zero elements of L
#'@param pi_theta Proportion of non-zero elements of theta
#'@param V_theta_prop Proportion of heritability contributed by theta. Length M
#'@param h2 Total heritability. Length M
#'@param R covariance of environmental effects. M times M
#'@export
sim_ind_data <- function(nind, nvar, F_mat, pi_L, pi_theta, V_theta_prop, h2, R){
  nfactor <- ncol(F_mat)
  ntrait <- nrow(F_mat)
  if(any(rowSums(F_mat==0) == nfactor)){stop("All traits must have a factor component")}
  stopifnot(length(V_theta_prop)==ntrait)
  stopifnot(length(h2)==ntrait)

  #Generate B_f = LF^T
  #Make it so that variance of loadings is 1
  sigma_L <- sqrt(1/pi_L)
  load_dist <- ashr::normalmix(pi=c(1-pi_L, pi_L), mean=rep(0, 2), sd=c(0, sigma_L))
  L_mat <- replicate(n=nfactor,causeSims::rnormalmix(nvar, load_dist))
  B_f <- L_mat%*%t(F_mat)
  #Variance from B_f (genos are standardizes)
  V_f <- colSums(B_f^2)

  #Generate theta
  V_theta <- V_theta_prop*V_f/(1-V_theta_prop)
  sigma_theta <- sqrt(V_theta/(pi_theta*nvar))
  theta <- map_dfc(sigma_theta, function(st){
    td <- ashr::normalmix(pi=c(1-pi_theta, pi_theta), mean=rep(0, 2), sd=c(0, st))
    causeSims::rnormalmix(nvar, td)
  }) %>% as.matrix()

  #Total genetic variance
  V_gen <- colSums((B_f + theta)^2)

  #Total environmental variance for each triat
  V_E <- (1-h2)*V_gen/h2
  #Covariance for environmental effects
  V_E_mat <- diag(sqrt(V_E)) %*% R %*%  diag(sqrt(V_E))

  #generate genotypes, scaled and centered
  maf <- rbeta(n=nvar, 1, 5)
  maf <- pmin(maf + 0.01, 0.99)
  G <- replicate(n = nind, rbinom(nvar, size=2, prob = maf)) %>% t()
  mu <- apply(G, 2, mean)
  s <- apply(G, 2, sd)
  G <- t((t(G)-mu)/s)

  #generate error
  E <- MASS::mvrnorm(n=nind, mu=rep(0, ntrait), Sigma=V_E_mat)

  M <- G%*%(B_f + theta) + E
  m <- apply(M, 2, mean)
  M <- t(t(M)-m)
  #Lazy summary statistics

  beta_hat <- crossprod(G, M)/(nind-1)
   #this is only approximate but it is good out to a few decimal places
  s <- sqrt(apply(M, 2, var)/(nind-1))
  se_hat <- matrix(rep(s, each=nvar), nrow=nvar, byrow=FALSE)

  ret <- list(L_mat = L_mat, F_mat = F_mat, theta = theta, G = G,
              E = E, M = M, sigma_theta = sigma_theta, sigma_L = sigma_L,
              beta_hat = beta_hat, se_hat = se_hat)
  #sumstats <- purrr::map(seq(ntrait), function(m){
  #                 b_m <- crossprod(G, M[,m])
  #                 x <- purrr::map(seq(nvar), function(j){
  #                      b_mj <- crossprod(G[,j], M[,m])/(nind-1)
  #                      s_mj <- sqrt( 1/(nind-2)*sum((M[,m]-b_mj*G[,j])^2)/(nind-1))
  #                      c(b_mj, s_mj)
  #                  })
  #                 b_m <- map(x, 1) %>% unlist()
  #                 s_m <- map(x, 2) %>% unlist()
  #})


}

if(FALSE){
library(sumstatFactors)
library(tidyverse)
library(flashier)
R1 <- matrix(0.7, nrow=4, ncol=4)
diag(R1) <- 1
R <- Matrix::bdiag(R1, R1, R1, R1, R1) %>% as.matrix()

R_ind <- diag(rep(1, 20))

F_mat <- readRDS("analysis_data/factors2.RDS")
F_mat[13,3] <- 0.3

 set.seed(1)
 #dat <- sim_ind_data(nind = 2000, nvar = 10000, F_mat = F_mat, pi_L = 0.1,
#                     pi_theta = 0.2, V_theta_prop = rep(0.3, 20), h2 = rep(0.1, 20), R = R)

 #dat <- sim_ind_data(nind = 2000, nvar = 10000, F_mat = F_mat, pi_L = 0.1,
 #                     pi_theta = 0.2, V_theta_prop = rep(0, 20), h2 = rep(0.1, 20), R = R)

 dat <- sim_ind_data(nind = 10000, nvar = 10000, F_mat = F_mat, pi_L = 0.1,
                     pi_theta = 0.2, V_theta_prop = rep(0, 20), h2 = rep(0.1, 20), R = R_ind)

 dat2 <- sim_bh3(true_L = dat$L_mat, true_F = dat$F_mat, true_Theta = dat$theta, S = dat$se_hat)

 fit_naive <- flash.init(dat2$beta_hat, S = dat2$se_hat, var.type = 2) %>%
              flash.add.greedy(., Kmax = 100, init.fn = init.fn.softImpute) %>%
              flash.backfit(.)

 ptrue <- plot_factors(F_mat, 1:20)
 pnaive <- plot_factors(fit_naive$loadings.pm[[2]], 1:20)


 #fit with z score should be similar
 fit_naive_z <- flash.init(data = with(dat2, beta_hat/se_hat), S = matrix(1, nrow=nrow(dat2$beta_hat), ncol=ncol(dat2$beta_hat)), var.type = 2) %>%
   flash.add.greedy(., Kmax = 100, init.fn = init.fn.softImpute) %>%
   flash.backfit(.)
 pnaive_z <- plot_factors(fit_naive_z$loadings.pm[[2]], 1:20)

 F_hat <- t(t(fit_naive_z$loadings.pm[[2]])*dat$se_hat[1,])
 F_mat_z <- F_mat/dat$se_hat[1,]

 #Rtilde <- cor(dat$M)
 R_eig <- eigen(R)
 U <- R_eig$vectors %*% diag(1/sqrt(R_eig$values))
 Z_tilde <- with(dat, beta_hat/se_hat) %*% U
 fit <- flash.init(data=Z_tilde, S = matrix(1, nrow=nrow(Z_tilde), ncol=ncol(Z_tilde)),  var.type=2) %>%
        flash.add.greedy(Kmax = 100, init.fn = init.fn.softImpute) %>%
        flash.backfit() %>%
        flash.nullcheck();
 F_hat <- R_eig$vectors %*% diag(sqrt(R_eig$values)) %*% fit$loadings.pm[[2]]
 pcorrected <- plot_factors(F_hat, 1:20)

 set.seed(2)
 nvar_null <- 1000
 nind <- nrow(dat$M)
 maf <- rbeta(n=nvar_null, 1, 5)
 maf <- pmin(maf + 0.01, 0.99)
 G_null <- replicate(n = nind, rbinom(nvar_null, size=2, prob = maf)) %>% t()
 mu <- apply(G_null, 2, mean)
 s <- apply(G_null, 2, sd)
 G_null <- t((t(G_null)-mu)/s)
 beta_hat_null <- crossprod(G_null, dat$M)/(nind-1)
 s <- sqrt(apply(dat$M, 2, var)/(nind-1))
 se_hat_null <- matrix(rep(s, each=nvar_null), nrow=nvar_null, byrow=FALSE)

 R_hat <- cor(beta_hat_null/se_hat_null)
 R_eig <- eigen(R_hat)
 U <- R_eig$vectors %*% diag(1/sqrt(R_eig$values))
 Z_tilde <- with(dat, beta_hat/se_hat) %*% U
 fit <- flash.init(data=Z_tilde, S = matrix(1, nrow=nrow(Z_tilde), ncol=ncol(Z_tilde)),  var.type=2) %>%
   flash.add.greedy(Kmax = 100, init.fn = init.fn.softImpute) %>%
   flash.backfit() %>%
   flash.nullcheck();
 F_hat <- R_eig$vectors %*% diag(sqrt(R_eig$values)) %*% fit$loadings.pm[[2]]
 pcorrected2 <- plot_factors(F_hat, 1:20)
}
