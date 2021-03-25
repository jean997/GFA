library(dplyr)
library(purrr)
library(sumstatFactors)
library(irlba)

args <- commandArgs(trailingOnly=TRUE)

out = args[1]
min_non_miss <- as.numeric(args[2])
p_thresh <- as.numeric(args[3])
k <- as.numeric(args[4])
flash_fit_file <- args[5]
snp_list_file <- args[6]
nb_files = args[-c(1:6)]


if(flash_fit_file == "NA"){

    # Read in data
    X <- map_dfr(nb_files, readRDS)
    snps <- readRDS(snp_list_file)
    X <- X %>% filter(snp %in% snps)
    ntrait <- X %>%
          select(ends_with(".est")) %>%
          ncol()
    nmiss <- X %>% 
         select(ends_with(".est")) %>%
         is.na() %>%
         rowSums()

    ix <- which(ntrait-nmiss >= min_non_miss)
    X <- X[ix,]

    Z_hat <- X %>% 
         select(ends_with(".est")) %>%
         as.matrix()

    snps <- X$snp
    nms <- names(X)[grep(".est$", names(X))]
}else{
    flash_fit <- readRDS(flash_fit_file)
    Z_hat <- flash_fit$fit$flash.fit$Y
    snps <- flash_fit$snps
    nms <- flash_fit$nms

}
if(p_thresh < 1){
    P <- 2*pnorm(-abs(Z_hat))
    Z_hat[P > p_thresh] <- 0
    n_small <- rowSums(Z_hat != 0)
    ix <- which(n_small > 0)
    Z_hat <- Z_hat[ix,]
    snps <- snps[ix]
}


f <-  irlba(Z_hat, nv = k)

f$snps <- snps
f$names <- nms


saveRDS(f, file=out)

