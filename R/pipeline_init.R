#'@title Set up Snakemake pipeline for summar statistic factors project
#'@export
pipeline_init <- function(){

  #Download code
  if(!dir.exists("R/")) system("mkdir R")
  files <- c("1_format_data.R","2_combine_data_vcf.R","ld_prune_one_chrom.R","ld_cat.R",
             "cause_params.R", "cause.R", "mrpresso.R", "lcv.R", "mr_package.R", "extract_results.R")
  for(f in files){
    download.file(url=paste0("https://gitlab.com/jean_morrison/gwas-summay-statistic-factorization/-/raw/master/pipeline/R/", f),
                  destfile = paste0("R/", f))
  }
  download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/Snakefile"),
                destfile = "Snakefile")
  download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/config.yaml"),
                destfile = "config.yaml")
  download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/cluster.yaml"),
                destfile = "cluster.yaml")
  download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/run-snakemake.sh"),
                destfile = "run-snakemake.sh")
  download.file(url="https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/gwas_pairs.csv",
                destfile="gwas_pairs.csv")
}
