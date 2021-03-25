library(dplyr)
library(purrr)
library(readr)


args <- commandArgs(trailingOnly=TRUE)
out <- args[1]
estL <- args[-1]

dat <- map_dfr(estL, function(f){
           eL <- readRDS(f)
           nf <- ncol(eL$L_est)
           bhat <- data.frame(eL$L_est)
           names(bhat) <- paste0("beta_", seq(nf))
           shat <- data.frame(eL$L_est_se)
           names(shat) <- paste0("se_", seq(nf))
           dat <- with(eL, cbind(P, bhat, shat))
           return(dat)
           })

write_tsv(dat, path=out)
