library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(readr)
library(purrr)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
c <- args[1]
gwas_info_file <- args[2]
norm_type <- args[3]
out <- args[4]
nm_out <- args[5]

info <- read_csv(gwas_info_file)

normbeta <- map(seq(nrow(info)),   function(i){
                        f <- info$formatted_data_path[i]
                        n <- info$name[i]
                        ss <- info$pub_sample_size[i]
                        if(is.na(ss)) ss <- info$pub_cases[i] + info$pub_controls[i]
                        v <- query_chrompos_file(paste0(c, ":1-536870911"), f)
                        dat <- vcf_to_tibble(v) %>% 
                                rename(snp=rsid)
                        if(norm_type=="ss"){
                            if(any(!is.na(dat$SS))){
                                dat <- dat %>% mutate(norm_beta  = (ES/SE)/sqrt(SS), 
                                                      norm_ss = 1/sqrt(SS)) %>%
                                    select(seqnames, start, snp, REF, ALT, norm_beta, norm_ss) 
                            }else{
                                dat <- dat %>% mutate(norm_beta  = (ES/SE)/sqrt(ss), 
                                                      norm_ss = 1/sqrt(ss)) %>%
                                    select(seqnames, start, snp, REF, ALT, norm_beta, norm_ss) 
                            }
                        }else if(norm_type == "z"){
                                dat <- dat %>% mutate(norm_beta  = ES/SE, 
                                                      norm_ss = 1) %>%
                                    select(seqnames, start, snp, REF, ALT, norm_beta, norm_ss) 
                        }
                        names(dat)[6] <- paste0(n, ".est")    
                        names(dat)[7] <- paste0(n, ".se")
                        return(dat)
                 }) %>%
            reduce(full_join, by = c("seqnames", "start", "snp", "REF", "ALT"))

t <- table(normbeta$snp)
snps <- names(t[t==1])
normbeta <- filter(normbeta, snp %in% snps)
            
saveRDS(normbeta, file=out)


# Save table of how traits are missing each SNP for LD clumping
miss <- normbeta %>% 
        select(-ends_with(".se")) %>%
        is.na(.) %>%
        rowSums(.)

nmiss <- data.frame(snp = normbeta$snp, miss = miss)
saveRDS(nmiss, file=nm_out)

