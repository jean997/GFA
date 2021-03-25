library(dplyr)
library(purrr)
library(stringr)
library(sumstatFactors)

args <- commandArgs(trailingOnly=TRUE)

out = args[1]
R_est_file <- args[2]
kmax <- as.numeric(args[3])
type <- args[4]
norm_type <- args[5]
seed <- as.numeric(args[6])
max_miss <- as.numeric(args[7])
p_thresh <- as.numeric(args[8])
n_sample <- as.numeric(args[9])
nb_files = args[-c(1:9)]

set.seed(seed)


stopifnot(type == "ev"| type =="plain" | str_starts(type, "ff-"))
if(str_starts(type, "ff")){
    max_ev <- as.numeric(str_replace(type, "ff-", ""))
    type <- "ff"
}

# Read in data
X <- map_dfr(nb_files, readRDS)

ntrait <- X %>%
          select(ends_with(".est")) %>%
          ncol()
nmiss <- X %>%
         select(ends_with(".est")) %>%
         is.na() %>%
         rowSums()

ix <- which(nmiss <= max_miss)
X <- X[ix,]

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
    R <- readRDS(R_est_file)
    stopifnot(all(R$names %in% nms))
    z_order <- match(R$names, nms)
    S_hat <- S_hat[,z_order]
    B_hat <- B_hat[,z_order]
    new_names <- R$names
   if(norm_type=="ss"){
        n <- 1/(S_hat[1,]^2)
    }
}else{
    new_names <- nms
}

if(n_sample < length(snps)){
    ix <- sample(seq_along(snps), size=n_sample, replace=FALSE)
    snps <- snps[ix]
    B_hat <- B_hat[ix,]
    S_hat <- S_hat[ix,]
}else{
    cat("SNP set is less than ", n_sample, ".\n")
}

if(norm_type == "ss"){
    if(type == "plain"){
        f <- fit_ff_new(B_std = B_hat, N = n, kmax = kmax)
    }else if(type == "ff"){
        f <- fit_ff_new(B_std = B_hat, N = n, R = R$R, kmax = kmax)
    }else if(type == "ev"){
        f <- fit_ev(B_hat, S_hat, n, R$R, kmax=kmax, adjust=FALSE)
    }
}else if(norm_type == "z"){
    if(type == "plain"){
        f <- fit_ff_new(Z_hat = B_hat, kmax = kmax)
    }else if(type == "ff"){
        f <- fit_ff_new(Z_hat = B_hat,R = R$R, kmax = kmax)
    }else if(type == "ev"){
        f <- fit_ev(B_hat, S_hat, n, R$R, kmax=kmax, adjust=FALSE, kmax = kmax)
    }
}

f$snps <- snps
f$names <- new_names

saveRDS(f, file=out)

