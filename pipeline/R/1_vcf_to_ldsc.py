library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]
sample_size <- args[3]
snp_ref <- args[4]



dat <- readVcf(input_file) %>%
           vcf_to_tibble()

ref <- read_table2(snp_ref)

if(all(!is.na(dat$SS))){  
     dat$N <- dat$SS
}else{
     dat$N <- sample_size
}
dat <- dat %>%
       rename(A1 = ALT, A2 = REF, SNP = rsid) %>%
       mutate(Z = ES/SE) %>%
       select(SNP, N, Z, A1, A2)

dat <- inner_join(dat, ref, by=c("SNP")) %>%
       mutate(A1_flip.x = recode(A1.x, C = "G", T = "A", A = "T", G = "C"), 
              A2_flip.x = recode(A2.x, C = "G", T = "A", A = "T", G = "C"), 
              sgn = case_when(A1.x == A1.y & A2.x == A2.y ~ 1, 
                              A1.x == A2.y & A2.x == A1.y ~ -1, 
                              A1_flip.x == A1.y & A2_flip.x == A2.y ~ 1, 
                              A1_flip.x == A2.y & A2_flip.x == A1.y ~ -1, 
                              TRUE ~ 0), 
              Z = Z*sgn,
              A1 = A1.y, A2 = A2.y) %>%
         select(SNP, N, Z, A1, A2)

write_tsv(dat, file=output_file)


