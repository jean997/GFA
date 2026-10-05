# susie_ss_ebnm_fn.R
#
# Same idea as susie_ebnm_fn.R (a flashier `ebnm_fn` that fits the loadings
# mode with a SuSiE regression in LD-whitened coordinates), but the per-block
# IBSS is done by susieR::susie_ss() instead of hand-written code.
#
# Per LD block b, with whitened design X_b (r_b x n_b; X_b'X_b = R_b), flashier
# gives us x (collapsed residual, length sum r_b) and constant s, and the
# regression is   x_b = X_b ell_b + N(0, s^2 I).   We call
#
#   susie_ss(XtX = R_b, Xty = X_b'x_b, yty = x_b'x_b, n = r_b,
#            standardize = FALSE, residual_variance = s^2,
#            estimate_residual_variance = FALSE,
#            prior_variance_grid = v, mixture_weights = pi)
#
# i.e. susieR's fixed-mixture-of-normals prior  beta ~ sum_k pi_k N(0, v_k).
#
# What this file still does itself (susieR does not provide it):
#   * the point mass at zero: it is represented by a numerically-null mixture
#     component with variance v0 = spike_rel * (smallest shat^2), because the
#     grid must be strictly positive;
#   * the EM update of the mixture weights pi, shared across blocks (susieR
#     treats the weights as fixed input);
#   * the ELBO handed back to flashier, assembled from fit$alpha, fit$mu,
#     fit$mu2 and fit$KL (so it does not depend on susieR's own ELBO bookkeeping);
#   * second-moment bookkeeping for flashier: Et2 = Et^2 + V/N.
#
# Same requirements as before: s constant across entries (column-wise variance
# type in flash_init, no missing data), flashier's L is X %*% ell (whitened).
#
# UNTESTED against a real susieR install (see test_susie_ss_ebnm.R).  Things to
# verify there are marked with "CHECK".

susie_ss_ebnm_fn <- function(design,
                             L = 5,
                             max_iter = 50,
                             tol = 1e-6,
                             grid_mult = sqrt(2),
                             spike_rel = 1e-10,
                             warm_start = TRUE,
                             solver = NULL) {
  if (is.null(solver)) {
    if (!requireNamespace("susieR", quietly = TRUE))
      stop("susie_ss_ebnm_fn needs the susieR package (>= 0.16, fixed mixture prior).")
    solver <- getExportedValue("susieR", "susie_ss")
  }
  if (is.matrix(design)) design <- list(design)
  info <- .ss_design_info(design)

  function(x, s, g_init = NULL, fix_g = FALSE,
           output = c("posterior_mean", "posterior_second_moment",
                      "fitted_g", "log_likelihood")) {
    inf <- info
    if (length(x) != info$N) inf <- .ss_design_info(list(diag(length(x))))  # flashier's startup test
    .ss_ebnm_call(x, s, g_init, fix_g, output, inf, L = L, max_iter = max_iter,
                  tol = tol, grid_mult = grid_mult, spike_rel = spike_rel,
                  warm_start = warm_start, solver = solver)
  }
}

# ---------------------------------------------------------------------------

.ss_design_info <- function(design) {
  design <- lapply(design, as.matrix)
  nr <- vapply(design, nrow, 1L); nc <- vapply(design, ncol, 1L)
  G  <- lapply(design, crossprod)
  list(X = design, G = G, d = lapply(G, diag), nr = nr, nc = nc,
       row_end = cumsum(nr), row_start = cumsum(nr) - nr + 1L,
       N = sum(nr), B = length(design))
}

.ss_lse_rows <- function(M) {
  mx <- M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
  mx + log(rowSums(exp(M - mx)))
}

.ss_sd_grid <- function(xhat, shat, mult) {
  smin <- min(shat) / 10
  v <- xhat^2 - shat^2
  smax <- if (any(v > 0)) 2 * sqrt(max(v)) else 8 * smin
  smax <- max(smax, 8 * smin)
  c(0, smin * mult^(0:ceiling(log(smax / smin) / log(mult))))
}

# Responsibility-weighted component counts for the EM update of pi, computed
# from the fitted alpha / mu and the residual of each effect (does not rely on
# susieR internals such as lbf_grid).
.ss_em_counts <- function(G, d, xtx_b, s2, alpha, mu, w, v) {
  K <- length(v); Lb <- nrow(alpha)
  Bm <- alpha * mu
  GB <- t(G %*% t(Bm)); Gmb <- colSums(GB)
  valid <- d > 1e-12 * max(d, 1e-300)
  nv <- sum(valid); shat2 <- s2 / d[valid]
  counts <- numeric(K)
  if (nv == 0) return(counts)
  for (l in seq_len(Lb)) {
    bh <- ((xtx_b - (Gmb - GB[l, ])) / d)[valid]
    lbf <- matrix(0, nv, K)
    for (k in seq_len(K))
      lbf[, k] <- -0.5 * log(1 + v[k] / shat2) + 0.5 * bh^2 * v[k] / (shat2 * (v[k] + shat2))
    lw <- sweep(lbf, 2, log(w), "+")
    r  <- exp(lw - .ss_lse_rows(lw))
    counts <- counts + colSums(alpha[l, valid] * r)
  }
  counts
}

.ss_ebnm_call <- function(x, s, g_init, fix_g, output, inf, L, max_iter, tol,
                          grid_mult, spike_rel, warm_start, solver) {
  N <- inf$N
  if (any(!is.finite(s)) || any(s <= 0)) stop("s must be finite and positive.")
  s0 <- s[1]
  if (max(abs(s - s0)) > 1e-6 * s0)
    stop("s must be constant across entries (column-wise variance type, no missing data).")
  s2 <- s0^2

  rows <- lapply(seq_len(inf$B), function(b) inf$row_start[b]:inf$row_end[b])
  xb   <- lapply(rows, function(r) x[r])
  xtx  <- lapply(seq_len(inf$B), function(b) as.vector(crossprod(inf$X[[b]], xb[[b]])))

  # ---- prior: sd grid (sd[1] = 0 is the null/spike component) and weights ----
  xh0 <- unlist(lapply(seq_len(inf$B), function(b) xtx[[b]] / pmax(inf$d[[b]], 1e-12)))
  sh0 <- unlist(lapply(seq_len(inf$B), function(b) s0 / sqrt(pmax(inf$d[[b]], 1e-12))))
  cand <- .ss_sd_grid(xh0, sh0, grid_mult)

  has_state <- inherits(g_init, "susie_ss_ebnm_state")
  if (has_state) {
    sd_k <- g_init$sd; pi_k <- g_init$pi
    if (!fix_g && max(sd_k) < 0.5 * max(cand)) {          # loading scale drifted up
      sd_k <- cand
      pi_k <- c(pi_k[1], rep((1 - pi_k[1]) / (length(sd_k) - 1), length(sd_k) - 1))
    }
  } else {
    sd_k <- cand
    pi_k <- c(0.9, rep(0.1 / (length(sd_k) - 1), length(sd_k) - 1))
  }
  pi_use <- pmax(pi_k, 1e-12); pi_use <- pi_use / sum(pi_use)
  v0 <- spike_rel * s2 / max(unlist(inf$d))
  v_all <- c(v0, sd_k[-1]^2)

  # ---- per-block susie_ss fits ----------------------------------------------
  prev <- if (has_state && warm_start) g_init$fits else NULL
  fits <- vector("list", inf$B); KLtot <- 0; Vtot <- 0
  Et <- numeric(N); counts <- numeric(length(v_all))

  for (b in seq_len(inf$B)) {
    Lb <- min(L, inf$nc[b])
    args <- list(XtX = inf$G[[b]], Xty = xtx[[b]], yty = sum(xb[[b]]^2),
                 n = max(inf$nr[b], 2L), L = Lb, standardize = FALSE,
                 residual_variance = s2, estimate_residual_variance = FALSE,
                 prior_variance_grid = v_all, mixture_weights = pi_use,
                 max_iter = max_iter, tol = tol, check_input = FALSE,
                 coverage = NULL, min_abs_corr = NULL)               # CHECK: NULL skips credible sets
    fit <- NULL
    if (!is.null(prev) && !is.null(prev[[b]]) && nrow(prev[[b]]$alpha) == Lb) {
      fit <- tryCatch(do.call(solver, c(args, list(model_init = prev[[b]]))),
                      error = function(e) NULL)                      # CHECK: model_init semantics
    }
    if (is.null(fit)) fit <- do.call(solver, args)

    alpha <- fit$alpha; mu <- fit$mu; mu2 <- fit$mu2
    G <- inf$G[[b]]; d <- inf$d[[b]]
    Bm <- alpha * mu
    GB <- t(G %*% t(Bm))
    ell <- colSums(Bm)
    Vb  <- sum(rowSums(alpha * mu2 * rep(d, each = nrow(alpha)))) - sum(Bm * GB)
    Et[rows[[b]]] <- as.vector(inf$X[[b]] %*% ell)
    KLtot <- KLtot + sum(fit$KL); Vtot <- Vtot + Vb
    counts <- counts + .ss_em_counts(G, d, xtx[[b]], s2, alpha, mu, pi_use, v_all)
    fits[[b]] <- list(alpha = alpha, mu = mu, mu2 = mu2, V = fit$V)
  }

  Et2  <- Et^2 + Vtot / N
  elbo <- -0.5 * (N * log(2 * pi * s2) + (sum((x - Et)^2) + Vtot) / s2) - KLtot

  # ELBO above is for the prior (sd_k, pi_use) that the fits used.  The EM
  # update below is what gets stored for the next call.
  pi_out <- if (fix_g) pi_k else (counts + 1e-10) / sum(counts + 1e-10)
  state <- structure(list(sd = sd_k, pi = pi_out, fits = fits, L = L),
                     class = "susie_ss_ebnm_state")

  if (identical(output, "lfsr"))
    return(list(posterior = data.frame(lfsr = rep(NA_real_, N)), fitted_g = state))
  if (identical(output, "posterior_sampler")) {
    sampler <- function(nsamp) {
      out <- matrix(0, nsamp, N)
      for (b in seq_len(inf$B)) {
        f <- fits[[b]]; ell <- matrix(0, nsamp, inf$nc[b])
        for (m in seq_len(nrow(f$alpha))) {
          idx <- sample.int(inf$nc[b], nsamp, replace = TRUE, prob = f$alpha[m, ])
          sdv <- sqrt(pmax(f$mu2[m, idx] - f$mu[m, idx]^2, 0))
          ell[cbind(seq_len(nsamp), idx)] <- ell[cbind(seq_len(nsamp), idx)] +
            rnorm(nsamp, f$mu[m, idx], sdv)
        }
        out[, rows[[b]]] <- ell %*% t(inf$X[[b]])
      }
      out
    }
    return(list(posterior_sampler = sampler, fitted_g = state))
  }

  list(posterior = data.frame(mean = Et, second_moment = Et2),
       fitted_g = state, log_likelihood = elbo, coefficients = Et)
}

# Per-variant PIPs and posterior means (original coordinates) per block.
susie_ss_state_summary <- function(state) {
  stopifnot(inherits(state, "susie_ss_ebnm_state"))
  lapply(state$fits, function(f)
    data.frame(pip = 1 - apply(1 - f$alpha, 2, prod), ell_pm = colSums(f$alpha * f$mu)))
}
