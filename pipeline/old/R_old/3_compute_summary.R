library(dplyr)
library(purrr)

args <- commandArgs(trailingOnly=TRUE)
normbeta <- readRDS(args[1])
keep_snps <- readRDS(args[2])
p_thresh <- as.numeric(args[3])
out_summ <- args[4]
out_data <- args[5]


nms <- names(normbeta)[grep(".est$", names(normbeta))]
n <- length(nms)

X <- normbeta %>%
      filter(snp %in% keep_snps$snp) 

#X <- right_join(normbeta, keep_snps, by="snp")
saveRDS(X, file=out_data)

B <-  X %>%
      #select(-snp, -seqnames, -start, -REF, -ALT) %>%
      select(ends_with(".est")) %>%
      as.matrix()

S <-  X %>%
      #select(-snp, -seqnames, -start, -REF, -ALT) %>%
      select(ends_with(".se")) %>%
      as.matrix()

pvals <- 2*pnorm(-abs(B/S))

dot_prds <- expand.grid(n1 = seq(n), n2 = seq(n)) %>%
            filter(n1 <= n2) 

prod <- apply(dot_prds, 1, function(x){
                        ix <- which(pvals[,x[1]] > p_thresh & 
                                    pvals[,x[2]] > p_thresh)
                        p <- B[ix,x[1]]*B[ix,x[2]]
                        s <- sum(p, na.rm=T)
                        nn <- sum(!is.na(p))
                        n1s <- sum(B[ix,x[1]][!is.na(p)])
                        n2s <- sum(B[ix,x[2]][!is.na(p)])
                        return(list(s = s, nn=nn, n1s=n1s, n2s=n2s))
                    })
dot_prds$prod <- map(prod, 1) %>% unlist()
dot_prds$n <- map(prod, 2) %>% unlist()
dot_prds$n1s <- map(prod, 3) %>% unlist()
dot_prds$n2s <- map(prod, 4) %>% unlist()
dot_prds$n1 <- nms[dot_prds$n1] %>% unlist()
dot_prds$n2 <- nms[dot_prds$n2] %>% unlist()

saveRDS(dot_prds, file=out_summ)




