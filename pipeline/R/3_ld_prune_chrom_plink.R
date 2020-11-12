library(tidyverse)
library(ieugwasr)

args <- commandArgs(trailingOnly=TRUE)
nmiss <- readRDS(args[1])
chrom <- as.numeric(args[2])
r2_thresh <- as.numeric(args[3])
ref_path  <- args[4]
out_nm = args[5]
#out_summ = args[7]

M <- max(nmiss$miss) + 1

nmiss <- nmiss %>% 
         mutate(pval = miss/M) %>%
         rename(rsid = snp)

nm_clump <- ld_clump(dat = nmiss, 
                     clump_r2 = r2_thresh, 
                     clump_p = 1, 
                     plink_bin = genetics.binaRies::get_plink_binary(), 
                     bfile = ref_path)

nm_clump <- rename(nm_clump, snp=rsid)

saveRDS(nm_clump, file=out_nm)

