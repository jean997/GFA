#'@title Convert GWAS summary statistics to a standard format
#'@description Format GWAS summary statistics for CAUSE
#'@param X data.frame
#'@param snp Column name containing SNP ID
#'@param beta_hat Column name containing effect estimate
#'@param se Column name containing standard error of beta_hat
#'@param A1 Column name containing effect allele
#'@param A2 Column name containing other allele
#'@param chrom Chromosome column (optional)
#'@param pos Position column (optional)
#'@param p_value p-value column (optional)
#'@param AF  Allele Frequency column (optional)
#'@param sample_size  Sample size column (optional) or an integer
#'@param compute_pval Logical, compute the p-value using a normal approximation if missing? Defaults to TRUE.
#'@param output_file File to write out formatted data. If missing formatted data will be returned.
#'@details This function will try to merge data sets X1 and X2 on the specified columns. Where
#'necessary, it will flip the sign of effects so that the effect allele is the same in both
#'data sets. It will remove variants with ambiguous alleles or where the alleles (G/C or A/T) or
#'with alleles that do not match between data sets (e.g A/G in one data set and A/C in the other).
#'It will not remove variants that are simply strand flipped between the two data sets (e. g. A/C in one data set, T/G in the other).
#'@return A data frame with columns chrom, pos, snp, A1, A2, beta_hat, se, p_value, and sample_size with all SNPs
#'aligned so that A is the effect allele. This is ready to be used with gwas_merge with formatted = TRUE.
#'@export
gwas_format <- function(X, snp, beta_hat, se, A1, A2,
                        chrom, pos, p_value, AF,
                        sample_size, output_file, compute_pval = TRUE){

  if(missing(snp) | missing(beta_hat) | missing(se) | missing(A1) | missing(A2)){
    stop("snp, beta_hat, se, A1, and A2 are required.\n")
  }
  if(missing(chrom)){
    X <- mutate(X, chrom = NA)
    chrom <- "chrom"
  }else if(is.na(chrom)){
    X <- mutate(X, chrom = NA)
    chrom <- "chrom"
  }
  if(missing(pos)){
    X <- mutate(X, pos = NA)
    pos <- "pos"
  }else if(is.na(pos)){
    X <- mutate(X, pos = NA)
    pos <- "pos"
  }
  if(missing(p_value)){
    X <- mutate(X, p_value = NA)
    p_value <- "p_value"
    p_val_missing <- TRUE
  }else if(is.na(p_value)){
    X <- mutate(X, p_value = NA)
    p_value <- "p_value"
    p_val_missing <- TRUE
  }else{
    p_val_missing <- FALSE
  }
  if(missing(AF)){
    X <- mutate(X, AF = NA)
    AF <- "AF"
  }else if(is.na(AF)){
    X <- mutate(X, AF = NA)
    AF <- "AF"
  }
  if(missing(sample_size)){
    X <- mutate(X, sample_size = NA)
    sample_size <- "sample_size"
  }else if(is.na(sample_size)){
    X <- mutate(X, sample_size = NA)
    sample_size <- "sample_size"
  }else if(is.numeric(sample_size)){
    X <- mutate(X, sample_size = sample_size)
    sample_size <- "sample_size"
  }
  keep_cols <- c(chrom, pos, snp, A1, A2, beta_hat, se, p_value, AF, sample_size)
  X <- X %>%
      select(keep_cols)%>%
      rename(snp = snp,
            beta_hat =beta_hat,
            se = se,
            A1 = A1,
            A2 = A2,
            chrom = chrom,
            pos = pos,
            p_value = p_value,
            AF = AF,
            sample_size = sample_size) %>%
      mutate(A1 = toupper(A1),
             A2 = toupper(A2))

  if(p_val_missing & compute_pval){
    X <- X %>% mutate(p_value = 2*pnorm(-abs(beta_hat/se)))
  }

  cat("There are ", nrow(X), " variants.\n")

  #Duplicated variants
  dup_vars <- X$snp[which(duplicated(X$snp))]
  X <- X %>% filter(!snp %in% dup_vars)
  cat("Removing ", length(dup_vars), " duplicated variants leaving ", nrow(X), "variants.\n")

  #Illegal alleles
  illegal_vars <- X %>%
                  filter((!A1 %in% c("A", "C", "T", "G") | !A2 %in% c("A", "C", "T", "G") )) %>%
                  select(snp)
  if(length(illegal_vars) > 0){
    X <- X %>% filter(!snp %in% illegal_vars$snp)
    cat("Removing ", length(illegal_vars), " variants with illegal alleles leaving ", nrow(X), "variants.\n")
  }else{
    cat("No variants have illegal alleles.\n")
  }

  #Ambiguous alleles
  n <- nrow(X)
  X <- remove_ambiguous(X, upper = TRUE)
  cat("Removed ", n-nrow(X), " variants with ambiguous strand.\n")

  cat("Flipping strand and effect allele so A1 is always A\n")
  X <- align_beta(X, "beta_hat", TRUE)


  X <- X %>% select(chrom, pos, snp, A1, A2, beta_hat, se, p_value, AF, sample_size)

  if(!missing(output_file)){
    cat("Writing out ", nrow(X), " variants to file.\n")
    write_tsv(X, path = output_file)
    return(NULL)
  }
  cat("Returning ", nrow(X), " variants.\n")
  return(X)

}

#'@export
read_standard_format <- function(file, ...){
  dat <- read_tsv(file, col_types=list(chrom="c",pos="i", A1 = "c", A2 = "c",
                       beta_hat="d", se = "d", p_value ="d", sample_size="d"), ...)
  return(dat)
}

#Flip signs and strabds so that allele 1 is allways A
align_beta <- function(X, beta_hat_name, upper=TRUE){
  flp = c("A" = "T", "G" = "C", "T" = "A",
          "C" = "G", "a"  = "t", "t" = "a",
          "c" = "g", "g" = "c")
  if(upper){
    X <- X %>% mutate( flip_strand = A1 == "T" | A2 == "T")
  }else{
    X <- X %>% mutate( flip_strand = A1 == "t" | A2 == "t")
  }
  X <- X %>% mutate(A1flp = case_when(flip_strand ~ flp[A1],
                                      TRUE ~ A1),
                    A2flp = case_when(flip_strand ~ flp[A2],
                                      TRUE ~ A2),
                    temp = case_when(A1flp == "A" | A1flp == "a" ~ get(beta_hat_name),
                                     TRUE ~ -1*get(beta_hat_name))) %>%
    select(-A1, -A2) %>%
    mutate(A1 = case_when(A1flp == "A" | A1flp=="a" ~ A1flp,
                          TRUE ~ A2flp),
           A2 = case_when(A1flp == "A" | A1flp=="a" ~ A2flp,
                          TRUE ~ A1flp)) %>%
    select(-A1flp, -A2flp, -flip_strand)
  ix <- which(names(X)==beta_hat_name)
  X <- X[,-ix]
  ix <- which(names(X) == "temp")
  names(X)[ix] <- beta_hat_name
  return(X)
}

remove_ambiguous <- function(X, upper=TRUE){
  if(upper){
    X <- X %>% dplyr::filter(!(A1 == "G" & A2 == "C") &
                               !(A1 == "C" & A2 == "G") &
                               !(A1 == "A" & A2 == "T") &
                               !(A1 == "T" & A2 == "A"))
    return(X)
  }
  X <- X %>% filter(!(A1 == "g" & A2 == "c") &
                      !(A1 == "c" & A2 == "g") &
                      !(A1 == "a" & A2 == "t") &
                      !(A1 == "t" & A2 == "a"))
  return(X)
}

