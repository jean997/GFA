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
                        chrom, pos, p_value,
                        sample_size, allele_freq,
                        output_file, compute_pval = TRUE){

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
    X <- mutate(X, pos = NA_integer_)
    pos <- "pos"
  }else if(is.na(pos)){
    X <- mutate(X, pos = NA_integer_)
    pos <- "pos"
  }

  if(missing(p_value)){
    X <- mutate(X, p_value = NA_real_)
    p_value <- "p_value"
    p_val_missing <- TRUE
  }else if(is.na(p_value)){
    X <- mutate(X, p_value = NA_real_)
    p_value <- "p_value"
    p_val_missing <- TRUE
  }else{
    p_val_missing <- FALSE
  }

  if(missing(sample_size)){
    X <- mutate(X, sample_size = NA_real_)
    sample_size <- "sample_size"
  }else if(is.na(sample_size)){
    X <- mutate(X, sample_size = NA_real_)
    sample_size <- "sample_size"
  }else if(is.numeric(sample_size)){
    X <- mutate(X, sample_size = sample_size)
    sample_size <- "sample_size"
  }

  if(missing(allele_freq)){
    X <- mutate(X, af = NA_real_)
    allele_freq <- "af"
  }else if(is.na(allele_freq)){
    X <- mutate(X, af = NA_real_)
    allele_freq <- "af"
  }

  keep_cols <- c(chrom, pos, snp, A1, A2, beta_hat, se, p_value, sample_size, allele_freq)

  X <- X %>%
      select(all_of(keep_cols))%>%
      rename(snp = snp,
            beta_hat =beta_hat,
            se = se,
            A1 = A1,
            A2 = A2,
            chrom = chrom,
            pos = pos,
            p_value = p_value,
            sample_size = sample_size,
            allele_freq = allele_freq) %>%
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

  # make X into a data.table for I can do the new align_beta on it.  ideally this would be for everything
  setDT(X)
  
  cat("Flipping strand and effect allele so A1 is always A\n")
  align_beta(X)

  # data table syntax
  X <- X[, .(chrom, pos, snp, A1, A2, beta_hat, se, p_value, sample_size, allele_freq)]

  if(!missing(output_file)){
    cat("Writing out ", nrow(X), " variants to file.\n")
    # changed from path= to file=
    write_tsv(X, file = output_file)
    return(NULL)
  }
  cat("Returning ", nrow(X), " variants.\n")
  return(X)

}


# read_standard_format <- function(file, ...){
#   dat <- read_tsv(file, col_types=list(chrom="c",pos="i", A1 = "c", A2 = "c",
#                        beta_hat="d", se = "d", p_value ="d", sample_size="d"), ...)
#   return(dat)
# }

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

# flip signs and strands so that allele 1 is always A
# now modifies X in-place w/ data table for speed and memory savings
# believe beta_hat and af are always the names assigned in gwas_format.  but kept for consistency w/ old func
align_beta <- function(X, upper = TRUE,
                          beta_col = "beta_hat",
                          af_col   = "af") {

  stopifnot(is.data.table(X))
  stopifnot(all(c("A1", "A2") %in% names(X)))
  stopifnot(beta_col %in% names(X))

  flp <- c("A"="T","G"="C","T"="A","C"="G",
           "a"="t","t"="a","c"="g","g"="c")

  af_present <- af_col %in% names(X)
  if (!af_present) {
    X[, (af_col) := NA_real_]  # create af col temporarily so code can be uniform
  }

  # flip strands if we have Ts to get As
  X[, flip_strands_flag :=
       if (upper) (A1 == "T" | A2 == "T") else (A1 == "t" | A2 == "t")]
  X[, `:=`(
    flipped_A1 = fifelse(flip_strands_flag, flp[A1], A1),
    flipped_A2 = fifelse(flip_strands_flag, flp[A2], A2)
  )]

  # flag that is true if the A1 is A (we want)
  X[, A1_A_flag := flipped_A1 %chin% c("A","a")]

  # swap A1 and A2 if A1 is not A
  X[, `:=`(
    A1 = fifelse(A1_A_flag, flipped_A1, flipped_A2),
    A2 = fifelse(A1_A_flag, flipped_A2, flipped_A1)
  )]

  # flip beta hats if the A1 is not A
  # sd means subset of data (select just beta column)
  X[, (beta_col) := {
    b <- .SD[[1L]]
    if (!is.numeric(b)) b <- as.numeric(b)
    fifelse(A1_A_flag, b, -b)
  }, .SDcols = beta_col]

  # flip af if the A1 is not A
  X[, (af_col) := {
    p <- .SD[[1L]]
    if (!is.numeric(p)) p <- as.numeric(p)
    fifelse(A1_A_flag, p, 1 - p)
  }, .SDcols = af_col]

  # delete columns we don't need anymore
  X[, c("flip_strands_flag", "flipped_A1", "flipped_A2", "A1_A_flag") := NULL]
  if (!af_present) X[, (af_col) := NULL]

  # since these are in-place mods, we can call func w/o assignment to new var.  but nobody wants to see the whole table
  return(invisible(X))
}
