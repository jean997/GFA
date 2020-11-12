library(dplyr)
library(purrr)
library(LaplacesDemon)

args <- commandArgs(trailingOnly=TRUE)
out = args[1]
files = args[-1]

x <- map(files, readRDS)


df <- map(x, function(y){mutate(y, prod = unlist(prod),
                                   ss = unlist(n),
                                   n1s = unlist(n1s),
                                   n2s  = unlist(n2s))}) %>%
        reduce(full_join, by=c("n1", "n2")) %>%
        mutate(prod = rowSums(select(., starts_with("prod"))),
               n = rowSums(select(., starts_with("ss"))),
               n1s = rowSums(select(., starts_with("n1s"))),
               n2s = rowSums(select(., starts_with("n2s")))) %>%
        select(n1, n2, prod, n, n1s, n2s) %>%
        mutate(cov = 1/(n-1)*(prod - n1s*n2s/n))

# we want a symmetric matrix so we need to replicate some rows of df
df_copy <- filter(df, n1 != n2) %>%
           rename(n1c = n2, n2c = n1) %>%
           rename(n1 = n1c, n2 = n2c)

cov_mat <- bind_rows(df, df_copy)  %>%
           select(n1, n2, cov) %>%
           reshape2::dcast(n1 ~ n2)

nms <- cov_mat$n1
R <- as.matrix(cov_mat[,-1]) %>%
     cov2cor() %>%
     as.symmetric.matrix() # Somehow cov2cor generates a matrix that is slightly non-symmetric

eS <- eigen(R)

ret <- list(R = R, names = nms, eS = eS)

saveRDS(ret, file=out)


