library(GFA)
library(flashier)

fixed_factors <- function(fit) {
  which(lengths(fit$flash_fit$fix.idx) > 0L)
}

# This signal causes unrestricted nullchecking to remove the fixed factor for R.
set.seed(401)
F <- matrix(c(1, 0.7, -0.6, 0.2), 4, 1)
Y <- matrix(rnorm(800), 800, 1) %*% t(3 * F) + matrix(rnorm(3200), 800, 4)
R <- matrix(0.2, 4, 4)
diag(R) <- 1
params <- gfa_default_parameters()
params$kmax <- 2

set.seed(10)
raw <- gfa_fit(Y, N = rep(100, 4), R = R, params = params, no_wrapup = TRUE)
stopifnot(is.null(raw$fit$flash_fit$maxiter.reached))
stopifnot(length(fixed_factors(raw$fit)) == 1L)

check_covariance <- function(fit) {
  fixed <- fixed_factors(fit)
  stopifnot(length(fixed) == 1L)
  W <- fit$F_pm[, fixed, drop = FALSE]
  # The low-rank representation has no truncation error in this example.
  stopifnot(isTRUE(all.equal(tcrossprod(W) + diag(0.8, 4), R,
                            tolerance = 1e-10)))
}

# Both public functions must retain the specified error covariance.
wrapped <- gfa_wrapup(raw$fit, method = "fixed_factors", scale = raw$scale,
                      nullcheck = TRUE)
check_covariance(wrapped$fit)
set.seed(10)
fit <- gfa_fit(Y, N = rep(100, 4), R = R, params = params)
check_covariance(fit$fit)
stopifnot(ncol(fit$F_hat) == fit$fit$n_factors - 1L)

# Backfitting can mark fixed factors as zero.
# With remove=TRUE, flash_nullcheck removes zero factors outside kset.
zero_fit <- flash_init(matrix(0, 20, 3), S = 1, var_type = NULL)
normal <- flash_ebnm(prior_family = "normal", scale = 1)
init <- list(matrix(1, 20, 1), matrix(1, 3, 1))
zero_fit <- flash_factors_init(zero_fit, init = init, ebnm_fn = normal)
zero_fit <- flash_factors_fix(zero_fit, kset = 1, which_dim = "factors")
zero_fit <- flash_factors_init(zero_fit, init = init, ebnm_fn = normal)
zero_fit <- flash_backfit(zero_fit, verbose = 0)
stopifnot(all(zero_fit$flash_fit$is.zero))
cleaned <- GFA:::gfa_nullcheck(zero_fit, tol = -Inf)
stopifnot(cleaned$n_factors == 1L, identical(fixed_factors(cleaned), 1L))
stopifnot(isTRUE(all.equal(cleaned$F_pm, zero_fit$F_pm[, 1, drop = FALSE])))

# An all-fixed fit has no candidates for removal. Retain its error covariance.
check_covariance(GFA:::gfa_nullcheck(flash_factors_remove(raw$fit, kset = 1)))

# Without fixed factors, the helper must match flashier's factor removal.
unfixed <- flash_factors_unfix(raw$fit, kset = fixed_factors(raw$fit))
pruned <- GFA:::gfa_nullcheck(unfixed)
stopifnot(pruned$n_factors == raw$fit$n_factors - 1L)
fields <- c("n_factors", "L_pm", "F_pm", "L_psd", "F_psd", "residuals_sd", "elbo")
stopifnot(isTRUE(all.equal(pruned[fields], flash_nullcheck(unfixed)[fields])))
unfixed <- flash_factors_unfix(zero_fit, kset = 1)
stopifnot(isTRUE(all.equal(GFA:::gfa_nullcheck(unfixed)[fields],
                          flash_nullcheck(unfixed)[fields])))
empty <- flash_init(matrix(0, 20, 3), S = 1)
stopifnot(GFA:::gfa_nullcheck(empty)$n_factors == 0L)
