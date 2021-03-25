library(dplyr)
library(purrr)
library(sumstatFactors)

args <- commandArgs(trailingOnly=TRUE)

out = args[1]
R_est_file <- args[2]
no_missing <- as.logical(args[3])
kmax <- as.numeric(args[4])
type <- args[5]
norm_type <- args[6]
seed <- as.numeric(args[7])
min_non_miss <- as.numeric(args[8])
p_thresh <- as.numeric(args[9])
nb_files = args[-c(1:9)]

set.seed(seed)

# Read in data
X <- map_dfr(nb_files, readRDS)

ntrait <- X %>%
          select(ends_with(".est")) %>%
          ncol()
nmiss <- X %>% 
         select(ends_with(".est")) %>%
         is.na() %>%
         rowSums()
if(no_missing){
    ix <- which(nmiss == 0)
    X <- X[ix,]
}else{
    ix <- which(ntrait-nmiss >= min_non_miss)
    X <- X[ix,]
}

B_hat <- X %>% 
         select(ends_with(".est")) %>%
         as.matrix()

S_hat <- X %>% 
         select(ends_with(".se")) %>%
         as.matrix()
snps <- X$snp
if(p_thresh < 1){
    P <- 2*pnorm(-abs(B_hat/S_hat))
    n_small <- rowSums(P < p_thresh, na.rm=TRUE)
    ix <- which(n_small > 0)
    B_hat <- B_hat[ix,]
    S_hat <- S_hat[ix,]
    snps <- snps[ix]
}

nms <- names(X)[grep(".est$", names(X))]

if(type=="ev" | type == "ff"){
   if(norm_type=="ss"){

        # For ff and ev we need every SNP to have the same sample size within traits
        # That is hard to get exactly so we will go for it approximately
        # This is only for ss type correction
        SS <- 1/S_hat^2
        # We will try just taking sample sizes within +/- 100 of the mode
        N <- apply(SS, 2, function(x){
           #t <- table(x)
           #n_mode <- as.numeric(names(t[which.max(t)]))
           n_mode <- max(x, na.rm=TRUE)
           ix <- which(abs(x-n_mode) <=2000)
           nix <- which(abs(x-n_mode) >2000)
           l <- length(ix)
           n <- mean(x[ix])
           return(list(n=n, ix=ix,nix=nix, l=l))
        })
        nvar <- nrow(X)
        S_new <- sapply(seq_along(N), function(i){
                x <- S_hat[,i]
                x[N[[i]]$nix] <- NA
                return(x)
                }) 
        B_new <- sapply(seq_along(N), function(i){
                x <- B_hat[,i]
                x[N[[i]]$nix] <- NA
                return(x)
                }) 
        nmiss <- rowSums(is.na(S_new))
        if(no_missing){
            ix <- which(nmiss == 0)
        }else{
            ix <- which(ntrait-nmiss >= min_non_miss)
        }
        S_new <- S_new[ix,]
        B_new <- B_new[ix,]
        snps <- snps[ix]
        S_hat <- S_new
        B_hat <- B_new
        R <- readRDS(R_est_file)
        n <- map(N, 1) %>% unlist() %>% as.numeric()
        stopifnot(all(R$names %in% nms))
        z_order <- match(R$names, nms)
        S_hat <- S_hat[,z_order]
        B_hat <- B_hat[,z_order]
    }else if(norm_type=="z"){
        R <- readRDS(R_est_file)
        n <- rep(1, ntrait)
        stopifnot(all(R$names %in% nms))
        z_order <- match(R$names, nms)
        S_hat <- S_hat[,z_order]
        B_hat <- B_hat[,z_order]
    }
}


if(type=="plain"){
 f <- fit_plain(B_hat, S_hat, adjust = FALSE, kmax=kmax)
}else if(type=="ev"){
 f <- fit_ev(B_hat, S_hat, n, R$R, kmax=kmax, adjust=FALSE)
}else if(type == "ff"){
 f <- fit_ff(B_hat, S_hat, n, R$R, kmax=kmax, adjust=FALSE)
}
f$snps <- snps
f$names <- nms

B_resid <- B_hat - f$B_hat
T_est <- apply(B_resid, 2, var)

f$T_est <- T_est

saveRDS(f, file=out)

