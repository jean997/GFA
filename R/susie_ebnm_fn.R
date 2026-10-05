# susie_ebnm_fn.R
#
# An `ebnm_fn` for flashier that replaces the n independent normal-means
# problems for a loadings mode with a SuSiE-style sparse regression in
# LD-whitened coordinates.
#
# Model, per LD block b and factor k (flashier hands us x and s for this k):
#     x_b = X_b ell_b + N(0, s^2 I),    X_b = D_b^{1/2} U_b'  (r_b x n_b),
# where R_b = U_b D_b U_b' is the (regularized / truncated) LD matrix, so
# X_b' X_b = R_b.  flashier must be run on the whitened data
#     Ytilde_b = D_b^{-1/2} U_b' Y_b,
# so that flashier's "L" is  Ltilde = X ell  (NOT ell itself; see
# susie_state_summary() for the per-variant loadings).
#
# ell_b = sum_{m=1}^L gamma_m beta_m  (L single effects),
# beta_m ~ g = sum_k pi_k N(0, sd_k^2), with sd_1 = 0 (point mass at zero),
# gamma_m ~ Uniform over the variants of the block.  g is shared across blocks.
#
# What is returned to flashier
#   posterior$mean          X %*% E[ell]
#   posterior$second_moment (X E[ell])^2 + V / N       (V = tr(R Cov(ell)), summed
#                                                      over blocks; spread evenly)
#   log_likelihood          the SuSiE ELBO for this regression
# flashier then computes  KL := log_likelihood - normal.means.loglik(x, s, Et, Et2),
# and because sum(Et2 - 2 x Et + x^2) = E||x - X ell||^2 exactly, this equals
# -KL(q || prior), which is the correct ELBO contribution.
#
# Requirements / assumptions
#   * s must be constant across entries (use a column-wise variance type in
#     flash_init, e.g. var_type = 2, and no missing data / exclusions), because
#     only the SUM of the second moments is meaningful.
#   * Base R only.  Dense R_b = crossprod(X_b) is stored per block.
#
# Usage sketch
#   design  <- lapply(blocks, function(R) { e <- eigen(R, symmetric = TRUE); keep <- ...;
#                     t(e$vectors[, keep] * rep(sqrt(e$values[keep]), each = nrow(R))) })
#   fn      <- susie_ebnm_fn(design, L = 5)
#   fl      <- flash_init(Ytilde, var_type = 2)
#   fl      <- flash_factors_init(fl, init, ebnm_fn = list(fn, ebnm_point_normal))
#   fl      <- flash_backfit(fl, extrapolate = FALSE)

susie_ebnm_fn <- function(design,
                          L = 5,
                          skip_z = NULL,
                          n_pi_updates = 2,
                          max_sweeps = 20,
                          alpha_tol = 1e-4,
                          grid_mult = sqrt(2)) {
  if (is.matrix(design)) design <- list(design)
  info <- .susie_design_info(design)

  function(x, s, g_init = NULL, fix_g = FALSE,
           output = c("posterior_mean", "posterior_second_moment",
                      "fitted_g", "log_likelihood")) {
    # flashier's startup check (test.ebnm.fn) calls us with a length-3 vector.
    # Fall back to an identity design in that case.
    inf <- info
    if (length(x) != info$N) inf <- .susie_design_info(list(diag(length(x))))

    .susie_ebnm_call(x, s, g_init, fix_g, output, inf, L = L, skip_z = skip_z,
                     n_pi_updates = n_pi_updates, max_sweeps = max_sweeps,
                     alpha_tol = alpha_tol, grid_mult = grid_mult)
  }
}

# ---------------------------------------------------------------------------

.susie_design_info <- function(design) {
  design <- lapply(design, as.matrix)
  nr <- vapply(design, nrow, 1L)
  nc <- vapply(design, ncol, 1L)
  G  <- lapply(design, crossprod)
  list(X = design, G = G, d = lapply(G, diag), nr = nr, nc = nc,
       row_end = cumsum(nr), row_start = cumsum(nr) - nr + 1L, N = sum(nr),
       B = length(design))
}

.lse_rows <- function(M) {
  mx <- M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
  mx + log(rowSums(exp(M - mx)))
}

.make_sd_grid <- function(xhat, shat, mult) {
  smin <- min(shat) / 10
  v <- xhat^2 - shat^2
  smax <- if (any(v > 0)) 2 * sqrt(max(v)) else 8 * smin
  smax <- max(smax, 8 * smin)
  n <- ceiling(log(smax / smin) / log(mult))
  c(0, smin * mult^(0:n))
}

# One single-effect regression with the shared mixture prior.
# Xtr: X'(x - X * other effects);  d = diag(X'X);  s2 = noise variance.
.ser_update <- function(Xtr, d, s2, pi_k, sd_k) {
  n <- length(Xtr)
  valid <- d > 1e-12 * max(d, 1e-300)
  alpha <- mu <- nu <- pnz <- numeric(n)
  counts <- numeric(length(pi_k))
  if (!any(valid)) {
    return(list(alpha = rep(1 / n, n), mu = mu, nu = nu, pnz = pnz, kl = 0, counts = counts))
  }
  xh  <- Xtr[valid] / d[valid]
  sh2 <- s2 / d[valid]
  sd2 <- sd_k^2
  V   <- outer(sh2, sd2, "+")                       # n_v x K
  ll  <- -0.5 * (log(2 * pi * V) + xh^2 / V)
  lw  <- sweep(ll, 2, log(pi_k), "+")
  lse <- .lse_rows(lw)                              # log m_g(xhat_i)
  r   <- exp(lw - lse)                              # component responsibilities
  pm  <- xh * sweep(1 / V, 2, sd2, "*")             # posterior mean | comp
  pv  <- sweep(sh2 / V, 2, sd2, "*")                # posterior var  | comp
  mu_v <- rowSums(r * pm)
  nu_v <- rowSums(r * (pv + pm^2))
  ll0  <- ll[, 1]                                   # sd_1 = 0: null likelihood
  lbf  <- lse - ll0
  nv   <- sum(valid)
  la   <- log(1 / nv) + lbf
  a    <- exp(la - max(la)); a <- a / sum(a)
  # KL(q(beta | gamma = i) || g), q exact posterior for N(xhat_i; beta, sh2_i)
  klmi <- -0.5 * log(2 * pi * sh2) - 0.5 * (xh^2 - 2 * xh * mu_v + nu_v) / sh2 - lse
  kl_alpha <- sum(ifelse(a > 0, a * log(a * nv), 0))
  alpha[valid] <- a; mu[valid] <- mu_v; nu[valid] <- nu_v
  pnz[valid] <- 1 - r[, 1]                        # P(beta != 0 | gamma = i)
  list(alpha = alpha, mu = mu, nu = nu, pnz = pnz,
       kl = kl_alpha + sum(a * klmi), counts = colSums(a * r))
}

# IBSS for one block until alpha stops changing.  Updates `st` in place (list).
.block_fit <- function(inf, b, st, xtx, s2, pi_k, sd_k, L, max_sweeps, alpha_tol) {
  G <- inf$G[[b]]; d <- inf$d[[b]]
  Bm  <- st$alpha * st$mu                                  # L x n
  GB  <- if (all(Bm == 0)) matrix(0, L, ncol(Bm)) else t(G %*% t(Bm))
  Gmb <- colSums(GB)
  KL  <- numeric(L); counts <- numeric(length(pi_k))
  for (sw in seq_len(max_sweeps)) {
    delta <- 0; counts[] <- 0
    for (m in seq_len(L)) {
      Xtr <- xtx - (Gmb - GB[m, ])
      u   <- .ser_update(Xtr, d, s2, pi_k, sd_k)
      delta <- max(delta, max(abs(u$alpha - st$alpha[m, ])))
      st$alpha[m, ] <- u$alpha; st$mu[m, ] <- u$mu; st$nu[m, ] <- u$nu
      st$pnz[m, ] <- u$pnz
      newGb <- as.vector(G %*% (u$alpha * u$mu))
      Gmb   <- Gmb - GB[m, ] + newGb
      GB[m, ] <- newGb
      KL[m] <- u$kl; counts <- counts + u$counts
    }
    if (delta < alpha_tol) break
  }
  Bm <- st$alpha * st$mu
  st$KL <- sum(KL); st$counts <- counts
  st$ell <- colSums(Bm)
  st$V <- sum(rowSums(st$alpha * st$nu * rep(d, each = L)) - rowSums(Bm * GB))
  st
}

.susie_ebnm_call <- function(x, s, g_init, fix_g, output, inf, L, skip_z,
                             n_pi_updates, max_sweeps, alpha_tol, grid_mult) {
  N <- inf$N
  if (any(!is.finite(s)) || any(s <= 0))
    stop("susie_ebnm_fn: s must be finite and positive (no exclusions).")
  s0 <- s[1]
  if (max(abs(s - s0)) > 1e-6 * s0)
    stop("susie_ebnm_fn: s must be constant across entries (use a column-wise ",
         "variance type, no missing data).")
  s2 <- s0^2

  xtx <- lapply(seq_len(inf$B), function(b)
    as.vector(crossprod(inf$X[[b]], x[inf$row_start[b]:inf$row_end[b]])))

  # --- warm start / prior ----------------------------------------------------
  has_state <- inherits(g_init, "susie_ebnm_state")
  blocks_ok <- has_state && !is.null(g_init$blocks) &&
    length(g_init$blocks) == inf$B &&
    all(vapply(seq_len(inf$B), function(b)
      identical(dim(g_init$blocks[[b]]$alpha), c(as.integer(L), inf$nc[b])), TRUE))

  xh0 <- unlist(lapply(seq_len(inf$B), function(b) xtx[[b]] / pmax(inf$d[[b]], 1e-12)))
  sh0 <- unlist(lapply(seq_len(inf$B), function(b) s0 / sqrt(pmax(inf$d[[b]], 1e-12))))
  cand <- .make_sd_grid(xh0, sh0, grid_mult)

  if (has_state && !is.null(g_init$sd)) {
    sd_k <- g_init$sd; pi_k <- g_init$pi
    if (!fix_g && max(sd_k) < 0.5 * max(cand)) {      # scale drifted: rebuild grid
      sd_k <- cand
      pi_k <- c(pi_k[1], rep((1 - pi_k[1]) / (length(sd_k) - 1), length(sd_k) - 1))
    }
  } else {
    sd_k <- cand
    pi_k <- c(0.9, rep(0.1 / (length(sd_k) - 1), length(sd_k) - 1))
  }
  stopifnot(sd_k[1] == 0)

  blocks <- vector("list", inf$B)
  for (b in seq_len(inf$B)) {
    if (blocks_ok) {
      blocks[[b]] <- g_init$blocks[[b]]
    } else {
      n <- inf$nc[b]
      blocks[[b]] <- list(alpha = matrix(1 / n, L, n), mu = matrix(0, L, n),
                          nu = matrix(0, L, n), pnz = matrix(0, L, n), KL = 0, counts = numeric(length(pi_k)),
                          ell = numeric(n), V = 0, skipped = FALSE)
    }
    blocks[[b]]$skipped <- FALSE
    if (!is.null(skip_z)) {
      z <- xtx[[b]] / (sqrt(pmax(inf$d[[b]], 1e-12)) * s0)
      if (max(abs(z)) < skip_z) blocks[[b]]$skipped <- TRUE
    }
  }

  fit_pass <- function(pi_k) {
    cnt <- numeric(length(pi_k))
    for (b in seq_len(inf$B)) {
      if (blocks[[b]]$skipped) {
        n <- inf$nc[b]
        blocks[[b]] <<- modifyList(blocks[[b]], list(
          alpha = matrix(1 / n, L, n), mu = matrix(0, L, n), nu = matrix(0, L, n),
          pnz = matrix(0, L, n), KL = -L * log(pi_k[1]), ell = numeric(n), V = 0))
        blocks[[b]]$counts <<- c(L, rep(0, length(pi_k) - 1))
      } else {
        blocks[[b]] <<- .block_fit(inf, b, blocks[[b]], xtx[[b]], s2, pi_k, sd_k,
                                   L, max_sweeps, alpha_tol)
      }
      cnt <- cnt + blocks[[b]]$counts
    }
    cnt
  }

  n_upd <- if (fix_g) 0 else n_pi_updates
  for (u in seq_len(n_upd)) {
    cnt  <- fit_pass(pi_k)
    pi_k <- (cnt + 1e-10) / sum(cnt + 1e-10)
  }
  fit_pass(pi_k)        # final pass at fixed pi: KL terms are valid for this g

  # --- assemble outputs --------------------------------------------------------
  Et <- numeric(N); KLtot <- 0; Vtot <- 0
  for (b in seq_len(inf$B)) {
    rows <- inf$row_start[b]:inf$row_end[b]
    Et[rows] <- as.vector(inf$X[[b]] %*% blocks[[b]]$ell)
    KLtot <- KLtot + blocks[[b]]$KL
    Vtot  <- Vtot + blocks[[b]]$V
  }
  Et2  <- Et^2 + Vtot / N
  elbo <- -0.5 * (N * log(2 * pi * s2) + (sum((x - Et)^2) + Vtot) / s2) - KLtot

  state <- structure(list(sd = sd_k, pi = pi_k, blocks = blocks, L = L),
                     class = "susie_ebnm_state")

  if (identical(output, "lfsr"))
    return(list(posterior = data.frame(lfsr = rep(NA_real_, N)), fitted_g = state))

  if (identical(output, "posterior_sampler")) {
    sampler <- function(nsamp) {
      out <- matrix(0, nsamp, N)
      for (b in seq_len(inf$B)) {
        bk <- blocks[[b]]; if (bk$skipped) next
        ell <- matrix(0, nsamp, inf$nc[b])
        for (m in seq_len(L)) {
          idx <- sample.int(inf$nc[b], nsamp, replace = TRUE, prob = bk$alpha[m, ])
          sdv <- sqrt(pmax(bk$nu[m, idx] - bk$mu[m, idx]^2, 0))  # moment-matched normal
          ell[cbind(seq_len(nsamp), idx)] <- ell[cbind(seq_len(nsamp), idx)] +
            rnorm(nsamp, bk$mu[m, idx], sdv)
        }
        out[, inf$row_start[b]:inf$row_end[b]] <- ell %*% t(inf$X[[b]])
      }
      out
    }
    return(list(posterior_sampler = sampler, fitted_g = state))
  }

  list(posterior = data.frame(mean = Et, second_moment = Et2),
       fitted_g = state, log_likelihood = elbo,
       coefficients = Et)   # lets coef() work in flashier's startup sign test
}

# Per-variant results (original, un-whitened coordinates) from a fitted state.
# Returns a list over blocks of data.frames with pip and posterior mean ell.
susie_state_summary <- function(state) {
  stopifnot(inherits(state, "susie_ebnm_state"))
  lapply(state$blocks, function(bk) {
    # PIP_i = 1 - prod_m (1 - alpha_mi * P(beta_m != 0 | gamma_m = i)): unused
    # effects (mass on the point-mass-at-zero component) then add ~nothing.
    data.frame(pip = 1 - apply(1 - bk$alpha * bk$pnz, 2, prod), ell_pm = bk$ell)
  })
}

# ---------------------------------------------------------------------------
# Helper: build the design (D^{1/2} U') and whitener (D^{-1/2} U') per LD block.
#   R_list    list of LD matrices, one per block
#   keep_frac eigenvalues below keep_frac * max(eigenvalue) are dropped
# Returns list(design = <list of r_b x n_b>, W = <list of r_b x n_b>).
# Whitened data for block b is  W[[b]] %*% Y[rows_of_block_b, ].
make_whitening <- function(R_list, keep_frac = 1e-3) {
  parts <- lapply(R_list, function(R) {
    e <- eigen(R, symmetric = TRUE)
    keep <- e$values > keep_frac * max(e$values)
    U <- e$vectors[, keep, drop = FALSE]; d <- e$values[keep]
    list(X = t(U * rep(sqrt(d), each = nrow(R))),
         W = t(U / rep(sqrt(d), each = nrow(R))))
  })
  list(design = lapply(parts, `[[`, "X"), W = lapply(parts, `[[`, "W"))
}
