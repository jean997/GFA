## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(sumstatFactors)
library(dplyr)
library(ggplot2)
library(viridis)

## ---- warning=FALSE-----------------------------------------------------------
set.seed(100)
myF <- matrix(c(0, -1, 1, 1, 2, 1), nrow = 3, byrow = TRUE)

myF <- apply(myF, 2, function(x){x/sum(x^2)})

p <- plot_factors(myF, names = c("T3", "T2", "T1"), factor_names = c("F1", "F2")) +
     scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2], name = "Normalized\nEffect") +
     xlab("Mediator") +
     ggtitle("True Mediator Network") +
     scale_x_discrete(position = "top")+
     theme(axis.text.x = element_text(angle = 0),
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank())
p

## -----------------------------------------------------------------------------
dat <- sim_sumstats_lf(F_mat = myF, N = 50000, J = 10000, h_2_trait = rep(0.5, 3),
                       omega = rep(0.8, 3), pi_L = rep(0.1,2), pi_theta = 0.1, 
                       overlap_prop =  0)
names(dat)

## ---- warning = FALSE---------------------------------------------------------

# Null SNPS which do not effect factors or have additional effects in theta
ix0 <- sample(which((rowSums(dat$L_mat == 0) == 2) & (rowSums(dat$theta == 0) == 3)), size = 8)
# SNPS which only effect trait 1 through theta and not through factors
ix1 <- sample(which((rowSums(dat$L_mat == 0) == 2) & (rowSums(dat$theta == 0) == 2) & abs(dat$theta[,1])> 0.02), size = 4)
# SNPS which only effect trait 2 through theta and not through factors
ix2 <- sample(which((rowSums(dat$L_mat == 0) == 2) & (rowSums(dat$theta == 0) == 2) & abs(dat$theta[,2])> 0.02), size = 5)
# SNPS which only effect trait 3 through theta and not through factors
ix3 <- sample(which((rowSums(dat$L_mat == 0) == 2) & (rowSums(dat$theta == 0) == 2) & abs(dat$theta[,3])> 0.02), size = 3)
# SNPS which only effect effect factor 1 only
ix4 <- sample(which((rowSums(dat$L_mat == 0) == 1) & abs(dat$L_mat[,1]) > 0.07 & (rowSums(dat$theta == 0) == 3)), size = 6)
# SNPS which effect factor 2 only 
ix5 <- sample(which((rowSums(dat$L_mat == 0) == 1) & abs(dat$L_mat[,2]) > 0.07 & (rowSums(dat$theta == 0) == 3)), size = 4)

X <- dat$beta_hat[c(ix0, ix1, ix2, ix3, ix4, ix5),][,c(3,2,1)]
se <- dat$se_beta_hat[c(ix0, ix1, ix2, ix3, ix4, ix5),][, c(3,2,1)]
plot_factors(X, factor_names = c("T1", "T2", "T3")) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2], name = "Beta hat") +
  xlab("Trait") +
  ylab("Variant") +
  ggtitle("Variant-Trait Associations") +
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 0),
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank())

## -----------------------------------------------------------------------------
R <- with(dat, est_R_pairwise(beta_hat, se_beta_hat, N = rep(100000, 3), z_scores = TRUE, p_val_thresh = 0.05))

## -----------------------------------------------------------------------------
plot_factors(R) + theme(axis.title =  element_blank())


## -----------------------------------------------------------------------------
eigen(R)$values

## -----------------------------------------------------------------------------
Z_hat = with(dat, beta_hat/se_beta_hat)
fit_plain <- fit_ff_prefit(Z_hat = Z_hat)

## -----------------------------------------------------------------------------
fit_withR <- fit_ff_prefit(Z_hat = Z_hat, R = R)

## -----------------------------------------------------------------------------
fit_withR$F_hat %>%
  plot_factors(., names = c("T3", "T2", "T1")) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2], name = "Normalized\nEffect") +
  xlab("Mediator") +
  ylab("Trait") +
  ggtitle("Estimated Mediator Network (GFA)") +
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 0),
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank())


## -----------------------------------------------------------------------------
fit_svd <- svd(Z_hat)

fit_svd$v %>%
  plot_factors(., names = c("T3", "T2", "T1")) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2], name = "Normalized\nEffect") +
  xlab("Mediator") +
  ylab("Trait") +
  ggtitle("Estimated Mediator Network (SVD)") +
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 0),
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank())


## -----------------------------------------------------------------------------
fit_dg <- sparse_svd(dat, z_score_thresh = 4)
fit_dg$v %>%
  plot_factors(., names = c("T3", "T2", "T1")) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2], name = "Normalized\nEffect") +
  xlab("Mediator") +
  ylab("Trait") +
  ggtitle("Estimated Mediator Network (DeGas)") +
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 0),
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank())

## -----------------------------------------------------------------------------
data("snpdata")
data("ld_evd_list")
head(snpdata)
length(ld_evd_list)

## -----------------------------------------------------------------------------
dim(snpdata)
#calculate the total number of SNPs with LD information. Verify it is the same as the number of SNPs in snpdata
sapply(ld_evd_list, function(x){length(x$values)}) %>% sum()

## -----------------------------------------------------------------------------
set.seed(200)
ntrait <- 50
nvar <- 200000
nfactor <- 12
N <- 50000
overlap_prop <- 1
pi_L <- pi_theta <- 5000/nvar

## -----------------------------------------------------------------------------
g_F <- function(n){runif(n, -1, 1)}
nz_factor <- c(pmin(rpois(3, 24)+1, ntrait), pmin(rpois(9, 4)+2, ntrait))

## -----------------------------------------------------------------------------
h_2_trait <- runif(n=ntrait, 0.05, 0.3)
omega <- runif(n=ntrait, 0.5, 1)
h_2_factor <- runif(n=nfactor, 0.5, 1)

## -----------------------------------------------------------------------------
Rblock <- matrix(0.3, 5, 5)
diag(Rblock) <- 1
R_E <- kronecker(diag(10), Rblock)

## -----------------------------------------------------------------------------
dat <-  sim_sumstats_lf(N=N, J=nvar, h_2_trait = h_2_trait,
            omega = omega, h_2_factor = h_2_factor,
            pi_L = rep(pi_L, nfactor), pi_theta = pi_theta,
            R_E = R_E, overlap_prop = overlap_prop,
            g_F = g_F, nz_factor = nz_factor, 
            R_LD = ld_evd_list, snp_info = snpdata) 

## -----------------------------------------------------------------------------
dat$F_mat %>%
  plot_factors(.) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2], name = "Normalized\nEffect") +
  xlab("Mediator") +
  ylab("Trait") +
  ggtitle("True Mediator Network") +
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 0),
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank())

## -----------------------------------------------------------------------------
dim(dat$snp_info)
head(dat$snp_info)

## -----------------------------------------------------------------------------
Z_hat <- with(dat, beta_hat/se_beta_hat)
ldscores <- dat$snp_info$ldscore_hm3
R <- R_ldsc(Z_hat, ldscores = ldscores, weights = 1/ldscores)

## -----------------------------------------------------------------------------
plot_factors(dat$R) + theme(axis.title =  element_blank()) + ggtitle("True R")
plot_factors(R) + theme(axis.title =  element_blank()) + ggtitle("Estimated R")
plot(dat$R, R, xlab = "True Correlation", ylab = "Estimated Correlation")
abline(0, 1)

## -----------------------------------------------------------------------------
ix <- which(dat$snp_info$keep_ld_prune_0.01);
Z_hat_indep <- Z_hat[ix,]
dim(Z_hat_indep)

## -----------------------------------------------------------------------------
fit <- fit_ff_prefit(Z_hat = Z_hat_indep, R = R)

## -----------------------------------------------------------------------------
dim(fit$F_hat)

## -----------------------------------------------------------------------------
fit$F_hat %>%
  plot_factors(.) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2], name = "Normalized\nEffect") +
  xlab("Mediator") +
  ylab("Trait") +
  ggtitle("Estimated Mediator Network (GFA)") +
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 0),
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank())

## -----------------------------------------------------------------------------
disc(fit$F_hat, dat$F_mat, lambda = 0)

## -----------------------------------------------------------------------------
fit_noR <- fit_ff_prefit(Z_hat = Z_hat_indep)

## -----------------------------------------------------------------------------
fit_noR$F_hat %>%
  plot_factors(.) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2], name = "Normalized\nEffect") +
  xlab("Mediator") +
  ylab("Trait") +
  ggtitle("Estimated Mediator Network (GFA with no R)") +
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 0),
           panel.background = element_rect(fill = "white"),
           axis.ticks = element_blank())

## -----------------------------------------------------------------------------
disc(fit_noR$F_hat, dat$F_mat, lambda = 0)

## -----------------------------------------------------------------------------
Z_hat_ht <- Z_hat_indep
Z_hat_ht[abs(Z_hat_ht) < 4] <- 0
fit_dg <- svd(Z_hat_ht)

disc(fit_dg$v, dat$F_mat, lambda = 0)

