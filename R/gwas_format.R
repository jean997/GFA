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
                        output_file, compute_pval = TRUE, return_og_snps = FALSE){

  # --- make data.table ---
  setDT(X)

  # --- check for missing inputs ---
  if(missing(snp) | missing(beta_hat) | missing(se) | missing(A1) | missing(A2)){
    stop("snp, beta_hat, se, A1, and A2 are required.\n")
  }
 
  if(missing(chrom) || is.na(chrom)){
    X[, chrom := NA]
    chrom <- "chrom"
  }
  if(missing(pos) || is.na(pos)){
    X[, pos := NA_integer_]
    pos <- "pos"
  }

  if(missing(p_value) || is.na(p_value)){
    X[, p_value := NA_real_]
    p_value <- "p_value"
    p_val_missing <- TRUE
  }else{
    p_val_missing <- FALSE
  }

  if(missing(sample_size) || is.na(sample_size)){
    X[, sample_size := NA_real_]
    sample_size <- "sample_size"
  }else if(is.numeric(sample_size)){
    X[, sample_size := NA_real_]
    sample_size <- "sample_size"
  }

  if(missing(allele_freq) || is.na(allele_freq)){
    X[, af := NA_real_]
    allele_freq <- "af"
  }
  
  # --- keep columns we want and rename ---
  old_cols <- c(chrom, pos, snp, A1, A2, beta_hat, se, p_value, sample_size, allele_freq)

  new_cols <- c(
    "chrom", "pos", "snp", "A1", "A2",
    "beta_hat", "se", "p_value", "sample_size", "allele_freq"
  )

  # checks
  stopifnot(all(old_cols %chin% names(X)))

  # drop columns not needed, by reference
  drop_cols <- setdiff(names(X), old_cols)
  if (length(drop_cols)) {
    X[, (drop_cols) := NULL]
  }

  # reorder columns
  setcolorder(X, old_cols)

  # rename columns, by reference
  setnames(X, old = old_cols, new = new_cols)

  # uppercase alleles, by reference
  X[, `:=`(
    A1 = toupper(A1),
    A2 = toupper(A2)
  )]

  cat("There are ", nrow(X), " variants.\n")

  # --- remove invalid snps ---
  # before filtering begins, label rows with indices, so we can join back to get all snps if they want
  # use ids instead of snp names bc may not be unique
  X[, row_id := .I]
  og_snp_index <- X[, .(row_id, snp)]
 
  #Duplicated variants
  # drop ALL instances of a variants which ever appears > 1x
  dup_vars <- X[duplicated(snp), unique(snp)]
  X <- X[!(snp %in% dup_vars)]
  cat("Removing", length(dup_vars), "duplicated variants leaving", nrow(X), "variants.\n")

  #Illegal alleles
  valid_alleles <- c("A", "C", "T", "G")
  illegal_vars <- X[
    !A1 %chin% valid_alleles | !A2 %chin% valid_alleles,
    snp
  ]
  if(length(illegal_vars)){
    print(head(X[snp %in% illegal_vars]))
    X <- X[!(snp %in% illegal_vars)]
    cat("Removing", length(illegal_vars), "variants with illegal alleles leaving", nrow(X), "variants.\n")
  }else{
    cat("No variants have illegal alleles.\n")
  }

  #Ambiguous alleles
  # when filtering rows, need to assign to var to keep result
  n <- nrow(X)
  X <- remove_ambiguous(X)
  cat("Removing", n-nrow(X), " variants with ambiguous strand leaving", nrow(X), "variants.\n")

  # --- compute pval ----
  # ask jean do we always want to compute even if not missing
  if(p_val_missing & compute_pval){
    X[, beta_hat := as.numeric(beta_hat)]
    X[, se       := as.numeric(se)]
    
    X[, p_value := fifelse(
      !is.na(beta_hat) & !is.na(se) & se != 0,
      2 * pnorm(-abs(beta_hat / se)),
      NA_real_
    )]
   cat("Computed p-value\n")
  }

  # --- harmonize alleles ---
  cat("Flipping strand and effect allele so A1 is always A\n")
  align_beta(X)

  # --- write out results ---
  if (return_og_snps){
    # label all snps that made it thru filtering as kept
    X[, pass_filt := TRUE]
    # make matrix with all snps given to function
    X_full <- X[
      og_snp_index,
      on = .(row_id, snp)
    ]
    # make NAs in kept from dropped rows into false
    X_full[is.na(pass_filt), pass_filt := FALSE]

    X <- X_full

    print(head(X))
  }
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

remove_ambiguous <- function(X) {
  stopifnot(data.table::is.data.table(X))
  stopifnot(all(c("A1", "A2") %chin% names(X)))

  ambig_pairs <- data.table(A1 = c("G", "C", "A", "T", "g", "c", "a", "t"),
                            A2 = c("C", "G", "T", "A", "c", "g", "t", "a"))

  idx_ambig <- X[ambig_pairs, on = .(A1, A2), which = TRUE, nomatch = NULL]

  print(head(X[idx_ambig]))

  if (!length(idx_ambig)) {
    return(invisible(X))
  } else {
    return(invisible(X[-idx_ambig]))
  }
}

# flip signs and strands so that allele 1 is always A
# now modifies X in-place w/ data table for speed and memory savings
# believe beta_hat and af are always the names assigned in gwas_format.  but kept for consistency w/ old func
align_beta <- function(X, beta_col = "beta_hat", af_col = "allele_freq") {

  # --- checks ---
  stopifnot(is.data.table(X))
  stopifnot(all(c("A1", "A2") %chin% names(X)))
  if (!is.character(X[["A1"]])) X[, A1 := as.character(A1)]
  if (!is.character(X[["A2"]])) X[, A2 := as.character(A2)]
  stopifnot(beta_col %chin% names(X))
  
  # --- setup ---
  flp <- c("A"="T","G"="C","T"="A","C"="G",
           "a"="t","t"="a","c"="g","g"="c")

  af_present <- af_col %chin% names(X)

  # make beta and af numeric if needed
  if (!is.numeric(X[[beta_col]])) {
    X[, (beta_col) := as.numeric(get(beta_col))]
  }

  if (af_present && !is.numeric(X[[af_col]])) {
    X[, (af_col) := as.numeric(get(af_col))]
  }

  # --- flipping ---
  # flip strands if we have Ts to get As
  idx_flip_strands <- X[, which(A1 %chin% c("T", "t") | A2 %chin% c("T", "t"))]

  if (length(idx_flip_strands)) {
    X[idx_flip_strands, `:=`(
      A1 = unname(flp[A1]),
      A2 = unname(flp[A2])
    )]
  }
  
  # we want A1 as A 
  idx_swap_alleles <- X[, which(!(A1 %chin% c("A", "a")))]
  
  if (length(idx_swap_alleles)) {
    X[idx_swap_alleles, `:=`(
      A1 = A2,
      A2 = A1
    )]

    # flip beta and af if A1 was not A
    X[idx_swap_alleles, (beta_col) := -get(beta_col)]
    if (af_present){
      X[idx_swap_alleles, (af_col)   := 1 - get(af_col)]
    }
  }

  # --- return ---
  # since these are in-place mods, we can call func w/o assignment to new var.  but nobody wants to see the whole table
  return(invisible(X))
}
