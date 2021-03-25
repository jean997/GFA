library(dplyr)
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

Z_hat <- X %>%
         select(ends_with(".est")) %>%
         as.matrix()

z_order <- match(R$names, nms)
Z_hat <- Z_hat[,z_order]

R$R <- LaplacesDemon::as.symmetric.matrix(R$R)

fake_fit <- list( fit = list(residuals.sd = rep(0, ncol(Z_hat))), F_hat = fit$v)


LL <- est_L_z(Z_hat = Z_hat, R = R$R, fit=fake_fit)
nfactor <- ncol(fake_fit$F_hat)
P <- with(LL, 2*pnorm(-abs(L_est/L_est_se)))
P <- cbind(X[, 1:5], P)
names(P)[-c(1:5)] <- paste0("p", 1:nfactor)
LL$P <- P

saveRDS(LL, out_file)
