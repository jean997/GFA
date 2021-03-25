library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
output <- args[1]
R_orig <- readRDS(args[2])
inp <- args[-c(1:2)]

res <- map_dfr(inp, function(f){
               #cat(f, "\n")
               r <- read_lines(f)
               n <- grep("Summary of Genetic Correlation Results", r)
               nms <- r[n+1] %>% str_squish() %>% str_split(" ") %>% unlist()
               d <- r[n+2] %>% str_squish() %>% str_split(" ") %>% unlist()
               d <- data.frame(val = d,  name = nms) %>%
                    spread(name, val)
        })

res$gcov_int <- as.numeric(res$gcov_int)
res$rg <- as.numeric(res$rg)
intercept_mat <- res %>%
                 mutate(n1 = str_replace(p1, "/project2/compbio/gwas_summary_statistics/standard_formats/ldsc/", "") %>%
                             str_replace(".tsv.gz", ""), 
                        n2 = str_replace(p2, "/project2/compbio/gwas_summary_statistics/standard_formats/ldsc/", "") %>%
                             str_replace(".tsv.gz", ""))  %>%
                 select(n1, n2, gcov_int) %>%
                 reshape2::dcast(n1~n2)

gc_mat <- res %>%
           mutate(n1 = str_replace(p1, "/project2/compbio/gwas_summary_statistics/standard_formats/ldsc/", "") %>%
                             str_replace(".tsv.gz", ""), 
                  n2 = str_replace(p2, "/project2/compbio/gwas_summary_statistics/standard_formats/ldsc/", "") %>%
                             str_replace(".tsv.gz", ""))  %>%
           select(n1, n2, rg) %>%
           reshape2::dcast(n1~n2)

saveRDS(res, file = output)
