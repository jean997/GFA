library(dplyr)
library(purrr)
library(sumstatFactors)

args <- commandArgs(trailingOnly=TRUE)
fit <- readRDS(args[1])
X <- readRDS(args[2])
R <- readRDS(args[3])
out_file <- args[4]

ntrait <- X %>%
          select(ends_with(".est")) %>%
          ncol()

nmiss <- X %>%
         select(ends_with(".est")) %>%
         is.na() %>%
         rowSums()

ix <- which(nmiss == 0)
X <- X[ix,]

nms <- names(X)[grep(".est$", names(X))]

B_hat <- X %>%
         select(ends_with(".est")) %>%
         as.matrix()

S_hat <- X %>%
         select(ends_with(".se")) %>%
         as.matrix()

z_order <- match(R$names, nms)
S_hat <- S_hat[,z_order]
B_hat <- B_hat[,z_order]

if(class(fit$fit$flash.fit$tau) == "numeric"){
    T_var <- 1/fit$fit$flash.fit$tau
}else{
    T_var <- 1/fit$fit$flash.fit$tau - fit$fit$flash.fit$given.S2
    T_var <- T_var[1,]
}

LL <- est_L(B_hat=B_hat, S_hat=S_hat, R = R$R, adjust=FALSE, tau = 1/T_var, fit=fit)
nfactor <- ncol(fit$F_hat)
P <- with(LL, 2*pnorm(-abs(L_est/L_est_se)))
P <- cbind(X[, 1:5], P)
names(P)[-c(1:5)] <- paste0("p", 1:nfactor)
LL$P <- P

saveRDS(LL, out_file)
