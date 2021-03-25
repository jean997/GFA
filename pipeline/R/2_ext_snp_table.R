library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(ieugwasr)
library(sumstatFactors)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)
c <- args[1]
orig_data_file <- args[2]
ext_phenos_file  <- args[3]
norm_type <- args[4]
out <- args[5]
nm_out <- args[6]

normbeta <- readRDS(orig_data_file)
trts <- read_table2(ext_phenos_file) # two cols triat and sample_size

snps <- normbeta$snp

phe <- map_df(1:ceiling(length(snps)/1000), function(i){
              k <- (i-1)*1000 + 1
              x <- associations(variants = snps[k:min(k + 999, length(snps))], id = trts$trait, proxies = 0)
              cat(k, " ")
              return(x)
              })
ref <- select(normbeta, snp, REF, ALT)
phe_aligned <- align_effects_ieu(phe, ref, ref_ea = "ALT", ref_nea = "REF")
phe_aligned$dup <- with(phe_aligned, duplicated(paste0(rsid, "-", id)))
phe_aligned <- phe_aligned %>% filter(!dup)

phe_aligned <- left_join(phe_aligned, trts, by = c("id" = "trait"))

#Merge old and new data
if(norm_type == "ss"){
    B_ext <- phe_aligned %>%
         filter(!dup) %>%
         mutate(trt = paste0(trait, "::", id), 
                est =  beta/(se*sqrt(sample_size))) %>%
         select(est, trt, rsid) %>%
         spread(trt, est)
    S_ext <- phe_aligned %>%
         filter(!dup) %>%
         mutate(trt = paste0(trait, "::", id), 
                se = case_when(!is.na(se) ~ 1/sqrt(sample_size), 
                               TRUE ~ se)) %>%
         select(se, trt, rsid) %>%
         spread(trt, se)
}else if(norm_type == "z"){
    B_ext <- phe_aligned %>%
         filter(!dup) %>%
         mutate(trt = paste0(trait, "::", id), 
                est =  beta/se) %>%
         select(est, trt, rsid) %>%
         spread(trt, est)
    S_ext <- phe_aligned %>%
         filter(!dup) %>%
         mutate(trt = paste0(trait, "::", id), 
                se =1) %>%
         select(se, trt, rsid) %>%
         spread(trt, se)
}
nb <- c("snp", paste0(names(B_ext)[-1], ".est"))
names(B_ext) <- nb
nms <- c(nms, nb[-1])
ns <- c("snp", paste0(names(S_ext)[-1], ".se"))
names(S_ext) <- ns

X <- full_join(normbeta, B_ext, by = "snp") %>%
     full_join(., S_ext, by = "snp") 

saveRDS(X, out)

# Save table of how traits are missing each SNP for LD clumping
miss <- X %>% 
        select(ends_with(".est")) %>%
        is.na(.) %>%
        rowSums(.)

nmiss <- data.frame(snp = normbeta$snp, miss = miss)
saveRDS(nmiss, file=nm_out)

