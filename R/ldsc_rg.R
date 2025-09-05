

################################################################################

#' Cross-Trait LD score regression
#'
#' @param ld_score Vector of LD scores.
#' @param ld_size Number of variants used to compute `ld_score`.
#' @param z1 Vector of z-scores for trait 1.
#' @param z2 Vector of z-scores for trait 2.
#' @param sample_size_1 Sample size of GWAS for trait 1.
#'   Possibly a vector, or just a single value.
#' @param sample_size_2 Sample size of GWAS for trait 2.
#'   Possibly a vector, or just a single value.
#' @param step1_chisq_max Threshold on `chi2` in step 1. Default is `30`.
#' @param chi2_thr2 Threshold on `chi2` in step 2. Default is `Inf` (none).
#' @param blocks Either a single number specifying the number of blocks,
#'   or a vector of integers specifying the block number of each `chi2` value.
#'   Default is `200` for `snp_ldsc()`, dividing into 200 blocks of approximately
#'   equal size. `NULL` can also be used to skip estimating standard errors,
#'   which is the default for `snp_ldsc2()`.
#' @param intercept You can constrain the intercept to some value (e.g. 0).
#'   Default is `NULL` (the intercept is estimated).
#'   Use a value of 0 if you are sure there is no overlap between GWAS samples.
#' @param intercept_h2_1 Intercept for heritability of trait 1 (default is NULL so the intercept is estimated).
#' @param intercept_h2_2 Intercept for heritability of trait 2 (default is NULL so the intercept is estimated).
#'
#' @return Vector of 4 values (or only the first 2 if `blocks = NULL`):
#'  - `[["int"]]`: LDSC regression intercept,
#'  - `[["int_se"]]`: SE of this intercept,
#'  - `[["h2"]]`: LDSC regression estimate of (SNP) heritability
#'  - `[["h2_se"]]`: SE of this heritability estimate.
#'
#'
#' @export
#' @keywords internal
ldsc_rg <- function(ld_score, ld_size, z1, z2, sample_size_1, sample_size_2,
                        blocks = NULL,
                        h2_1 = NULL,
                        h2_2 = NULL,
                        intercept = NULL,
                        intercept_h2_1 = NULL,
                        intercept_h2_2 = NULL,
                        step1_chisq_max = 30,
                        chi2_thr2 = Inf,
                        ncores = 1) {

  M <- length(z1)
  stopifnot(length(z2) == M)
  stopifnot(length(ld_score) == M)

  if (length(sample_size_1) == 1) {
    sample_size_1 <- rep(sample_size_1, M)
  } else {
    stopifnot(length(sample_size_1) == M)
  }
  if (length(sample_size_2) == 1) {
    sample_size_2 <- rep(sample_size_2, M)
  } else {
    stopifnot(length(sample_size_2) == M)
  }

  # First compute heritabilities
  # if(h2_se){
  #   if(is.null(blocks)) stop("If h2_se is desired, please provide blocks.\n")
  #   h2_blocks <- blocks
  # }else{
  #   h2_blocks <- NULL
  # }

  if(is.null(h2_1)){
    h2_1 <- snp_ldsc(ld_score, ld_size, z1^2,
                   sample_size_1,
                   blocks = NULL,
                   chi2_thr1 = step1_chisq_max,
                   chi2_thr2 = chi2_thr2)
  }
  pred_h2_1 <- h2_1[["int"]] + h2_1[["h2"]]*sample_size_1*ld_score/ld_size


  if(is.null(h2_2)){
    h2_2 <- snp_ldsc(ld_score, ld_size, z2^2,
                   sample_size_2,
                   blocks = NULL,
                   chi2_thr1 = step1_chisq_max,
                   chi2_thr2 = chi2_thr2)
  }
  pred_h2_2 <- h2_2[["int"]] + h2_2[["h2"]]*sample_size_2*ld_score/ld_size

  step1_index <- which(z1^2 < step1_chisq_max & z2^2 < step1_chisq_max)

  if(is.null(blocks)){
    result <- snp_ldsc(ld_score, ld_size, z1*z2,
                       sample_size = sqrt(sample_size_1*sample_size_2),
                       blocks = NULL, #blocks,
                       intercept= intercept,
                       chi2_thr2 = chi2_thr2,
                       step1_index = step1_index,
                       w0 = pred_h2_1*pred_h2_2,
                       type = "rg")

    if(h2_1[["h2"]] < 0 | h2_2[["h2"]] < 0){
      warning("Negative heritability estimates, genetic correlation will be missing.")
      gencorr <- NA
    }else{
      gencorr <- result[["h2"]]/sqrt(h2_1[["h2"]]*h2_2[["h2"]])
    }

    res <- c(int = result[["int"]],
             gencov = result[["h2"]],
             t1_int = h2_1[["int"]],
             t2_int = h2_2[["int"]],
             h2_1 = h2_1[["h2"]],
             h2_2 = h2_2[["h2"]],
             gencorr = gencorr)
  }else if(!is.null(blocks)){
      #### delete-a-group jackknife variance estimator ####
      if (length(blocks) == 1) {
        blocks <- sort(rep_len(seq_len(blocks), M))
      } else {
        stopifnot(length(blocks) == M)
      }
      ind_blocks <- split(seq_along(blocks), blocks)
      h_blocks <- M / lengths(ind_blocks)

      bigparallelr::register_parallel(ncores)

      delete_values <- foreach(ind_rm = c(list(NULL), ind_blocks), .combine = "cbind") %dopar% {
        keep <- which(!seq_along(ld_score) %in% ind_rm)
        ldsc_rg(ld_score = ld_score[keep],
                            ld_size,
                            z1 = z1[keep], z2 = z2[keep],
                            sample_size_1 = sample_size_1[keep],
                            sample_size_2 = sample_size_2[keep],
                            blocks = NULL,
                            h2_1 = NULL,
                            h2_2 = NULL,
                            intercept = NULL,
                            intercept_h2_1 = NULL,
                            intercept_h2_2 = NULL,
                            step1_chisq_max = step1_chisq_max,
                            chi2_thr2 = chi2_thr2,
                            ncores = 1)
      }
      estim <- delete_values[, 1]

      # https://doi.org/10.1023/A:1008800423698
      int_pseudovalues <- h_blocks * estim[1] - (h_blocks - 1) * delete_values[1, -1]
      gencov_pseudovalues  <- h_blocks * estim[2] - (h_blocks - 1) * delete_values[2, -1]
      h21_pseudovalues  <- h_blocks * estim[5] - (h_blocks - 1) * delete_values[5, -1]
      h22_pseudovalues  <- h_blocks * estim[6] - (h_blocks - 1) * delete_values[6, -1]
      gencorr_pseudovalues  <- h_blocks * estim[7] - (h_blocks - 1) * delete_values[7, -1]

      int_J <- sum(int_pseudovalues / h_blocks)
      gencov_J <- sum(gencov_pseudovalues / h_blocks)
      h21_J  <- sum( h21_pseudovalues / h_blocks)
      h22_J  <- sum( h22_pseudovalues / h_blocks)
      gencorr_J <- sum(gencorr_pseudovalues/ h_blocks)
      res <- c(int = int_J,
               int_se = sqrt(mean((int_pseudovalues - int_J)^2 / (h_blocks - 1))),
               gencov = gencov_J,
               gencov_se  = sqrt(mean(( gencov_pseudovalues -  gencov_J)^2 / (h_blocks - 1))),
               t1_int = h2_1[["int"]],
               t2_int = h2_2[["int"]],
               h2_1     = h21_J,
               h2_1_se  = sqrt(mean(( h21_pseudovalues -  h21_J)^2 / (h_blocks - 1))),
               h2_2     = h22_J,
               h2_2_se  = sqrt(mean(( h21_pseudovalues -  h22_J)^2 / (h_blocks - 1))),
               gencorr    = gencorr_J,
               gencorr_se  = sqrt(mean(( gencorr_pseudovalues -  gencorr_J)^2 / (h_blocks - 1))))
  }

  return(res)
}
