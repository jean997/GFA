################################################################################

# heteroscedasticity and overcounting weights
WEIGHTS_h2 <- function(pred, w_ld){
  1 / ((2*pred^2) *w_ld)
}

WEIGHTS_rg <- function(pred, w_ld, w0){
  od_w <- w0 + pred^2
  1 / (od_w*w_ld)
}

crossprod2 <- function(x, y) drop(base::crossprod(x, y))

# equivalent to stats::lm.wfit(cbind(1, x), y, w)
wlm <- function(x, y, w) {
  wx <- w * x
  W   <- sum(w)
  WX  <- sum(wx)
  WY  <- crossprod2(w,  y)
  WXX <- crossprod2(wx, x)
  WXY <- crossprod2(wx, y)
  alpha <- (WXX * WY - WX * WXY) / (W * WXX - WX^2)
  beta  <- (WXY * W  - WX * WY)  / (W * WXX - WX^2)
  list(intercept = alpha, slope = beta, pred = x * beta + alpha)
}

# equivalent to stats::lm.wfit(as.matrix(x), y, w)
wlm_no_int <- function(x, y, w) {
  wx <- w * x
  WXX <- crossprod2(wx, x)
  WXY <- crossprod2(wx, y)
  beta  <- WXY / WXX
  list(slope = beta, pred = x * beta)
}

################################################################################

#' LD score regression
#'
#' @param ld_score Vector of LD scores.
#' @param ld_size Number of variants used to compute `ld_score`.
#' @param chi2 Vector of chi-squared statistics.
#' @param sample_size Sample size of GWAS corresponding to chi-squared statistics.
#'   Possibly a vector, or just a single value.
#' @param chi2_thr1 Threshold on `chi2` in step 1. Default is `30`.
#'   This is equivalent to parameter `--two-step`.
#' @param chi2_thr2 Threshold on `chi2` in step 2. Default is `Inf` (none).
#' @param blocks Either a single number specifying the number of blocks,
#'   or a vector of integers specifying the block number of each `chi2` value.
#'   Default is `200` for `snp_ldsc()`, dividing into 200 blocks of approximately
#'   equal size. `NULL` can also be used to skip estimating standard errors,
#'   which is the default for `snp_ldsc2()`.
#' @param intercept You can constrain the intercept to some value (e.g. 1).
#'   Default is `NULL` in `snp_ldsc()` (the intercept is estimated)
#'   and is `1` in `snp_ldsc2()` (the intercept is fixed to 1).
#'   This is equivalent to parameter `--intercept-h2`.
#' @param type,w0,step1_index parameters used when called from snp_ldsc_rg
#'
#' @return Vector of 4 values (or only the first 2 if `blocks = NULL`):
#'  - `[["int"]]`: LDSC regression intercept,
#'  - `[["int_se"]]`: SE of this intercept,
#'  - `[["h2"]]`: LDSC regression estimate of (SNP) heritability (also see
#'    [coef_to_liab]),
#'  - `[["h2_se"]]`: SE of this heritability estimate.
#'
#'
#' @export
#' @keywords internal
snp_ldsc <- function(ld_score, ld_size, chi2, sample_size,
                     blocks = 200,
                     intercept = NULL,
                     chi2_thr1 = 30,
                     chi2_thr2 = Inf,
                     ncores = 1,
                     step1_index = NULL,
                     type = c("h2", "rg"),
                     w0= NULL){

  type <- match.arg(type, c("h2", "rg"))
  M <- length(chi2)
  if(type == "h2"){
    chi2 <- chi2 + 1e-8
    stopifnot(all(chi2 > 0))
  }else if(type == "rg"){
    stopifnot(length(w0) == M)
  }
  stopifnot(length(ld_score) == M)


  if (length(sample_size) == 1) {
    sample_size <- rep(sample_size, M)
  } else {
    stopifnot(length(sample_size) == M)
  }

  if (is.null(blocks)) {

    #### step 1 ####

    step1_int <- if (is.null(intercept)) {
      if(!is.null(step1_index)){
        ind_sub1 <- step1_index
      }else{
        ind_sub1 <- which(chi2 < chi2_thr1)
      }
      w_ld <- pmax(ld_score[ind_sub1], 1)
      x1 <- (ld_score / ld_size * sample_size)[ind_sub1]
      y1 <- chi2[ind_sub1]

      if(type == "h2"){
        pred0 <- y1
        wt_fun <- WEIGHTS_h2
      }else if(type == "rg"){
        pred0 <- 1/w0[ind_sub1]
        wt_fun <- function(pred, w_ld){ WEIGHTS_rg(pred, w_ld, w0 = w0[ind_sub1])}
      }
      for (i in 1:100) {
        pred <- wlm(x1, y1, w = wt_fun(pred0, w_ld))$pred
        if (max(abs(pred - pred0)) < 1e-6) break
        pred0 <- pred
      }
      wlm(x1, y1, w = wt_fun(pred0, w_ld))$intercept

    } else intercept

    #### step 2 ####

    ind_sub2 <- which(chi2 < chi2_thr2)
    w_ld <- pmax(ld_score[ind_sub2], 1)
    x <- (ld_score / ld_size * sample_size)[ind_sub2]
    y <- chi2[ind_sub2]
    yp <- y - step1_int

    if(type == "h2"){
      pred0 <- y
      wt_fun <- WEIGHTS_h2
    }else if(type == "rg"){
      pred0 <- 1/w0[ind_sub2]
      wt_fun <- function(pred, w_ld){ WEIGHTS_rg(pred, w_ld, w0 = w0[ind_sub2])}
    }
    for (i in 1:100) {
      pred <- step1_int + wlm_no_int(x, yp, w = wt_fun(pred0, w_ld))$pred
      if (max(abs(pred - pred0)) < 1e-6) break
      pred0 <- pred
    }
    step2_h2 <- wlm_no_int(x, yp, w = wt_fun(pred0, w_ld))$slope

    c(int = step1_int, h2 = step2_h2)

  } else {

    #### delete-a-group jackknife variance estimator ####

    if (length(blocks) == 1) {
      blocks <- sort(rep_len(seq_len(blocks), M))
    } else {
      stopifnot(length(blocks) == M)
    }
    ind_blocks <- split(seq_along(blocks), blocks)
    h_blocks <- M / lengths(ind_blocks)

    bigparallelr::register_parallel(ncores)

    s1i <- NULL
    delete_values <- foreach(ind_rm = c(list(NULL), ind_blocks), .combine = "cbind") %dopar% {
      keep <- which(!seq_along(ld_score) %in% ind_rm)
      if(!is.null(step1_index)) s1i <- which(keep %in% step1_index)
      snp_ldsc(ld_score[keep], ld_size, chi2[keep], sample_size[keep],
               NULL, intercept, chi2_thr1, chi2_thr2,
               1, s1i, type, w0[keep])
    }
    estim <- delete_values[, 1]

    # https://doi.org/10.1023/A:1008800423698
    int_pseudovalues <- h_blocks * estim[1] - (h_blocks - 1) * delete_values[1, -1]
    h2_pseudovalues  <- h_blocks * estim[2] - (h_blocks - 1) * delete_values[2, -1]

    int_J <- sum(int_pseudovalues / h_blocks)
    h2_J  <- sum( h2_pseudovalues / h_blocks)

    c(int    = int_J,
      int_se = sqrt(mean((int_pseudovalues - int_J)^2 / (h_blocks - 1))),
      h2     = h2_J,
      h2_se  = sqrt(mean(( h2_pseudovalues -  h2_J)^2 / (h_blocks - 1))))

  }
}
