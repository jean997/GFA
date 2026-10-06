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
# ell_b = sum_{m=1}^{L_b} gamma_m beta_m  (L_b = min(L, n_b) single effects, so a
# block never has more effects than variants; blocks with one variant and one
# whitened row use a vectorized closed-form path, see .singleton_update),
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
                          grid_mult = sqrt(2),
                          fast_singletons = TRUE) {
  if (is.matrix(design)) design <- list(design)
  info <- .susie_design_info(design, fast_singletons)

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

# Blocks with a single variant and a single (whitened) row are "singletons".
# They are handled together by a vectorized closed-form update (one effect on
# one variant, alpha = 1) instead of an IBSS loop; set fast_singletons = FALSE
# to send them through the general block path (used for testing).
.susie_design_info <- function(design, fast_singletons = TRUE) {
  design <- lapply(design, as.matrix)
  B  <- length(design)
  nr <- vapply(design, nrow, 1L)
  nc <- vapply(design, ncol, 1L)
  single <- if (fast_singletons) which(nr == 1L & nc == 1L) else integer(0)
  multi  <- setdiff(seq_len(B), single)
  G <- d <- vector("list", B)
  for (b in multi) { G[[b]] <- crossprod(design[[b]]); d[[b]] <- diag(G[[b]]) }
  sx <- vapply(design[single], function(M) M[1, 1], 0)       # numeric(0) if none
  row_end <- cumsum(nr)
  list(X = design, G = G, d = d, nr = nr, nc = nc,
       row_end = row_end, row_start = row_end - nr + 1L, N = sum(nr), B = B,
       single = single, multi = multi,
       sx = sx, sd2 = sx^2, srow = (row_end - nr + 1L)[single])
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

# Same calculation for MANY singleton blocks at once.  Each has one effect on
# one variant (alpha = 1, so the alpha-KL term is 0), so this is just the
# per-variant posterior / marginal likelihood / KL under g, vectorized over blocks.
# xh, sh2: marginal estimate and its variance for each (active) singleton.
.singleton_update <- function(xh, sh2, pi_k, sd_k) {
  sd2 <- sd_k^2
  V   <- outer(sh2, sd2, "+")
  ll  <- -0.5 * (log(2 * pi * V) + xh^2 / V)
  lw  <- sweep(ll, 2, log(pi_k), "+")
  lse <- .lse_rows(lw)
  r   <- exp(lw - lse)
  pm  <- xh * sweep(1 / V, 2, sd2, "*")
  pv  <- sweep(sh2 / V, 2, sd2, "*")
  mu  <- rowSums(r * pm)
  nu  <- rowSums(r * (pv + pm^2))
  kl  <- -0.5 * log(2 * pi * sh2) - 0.5 * (xh^2 - 2 * xh * mu + nu) / sh2 - lse
  list(mu = mu, nu = nu, pnz = 1 - r[, 1], kl = kl, counts = colSums(r))
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
  N <- inf$N; B <- inf$B; single <- inf$single; multi <- inf$multi
  if (any(!is.finite(s)) || any(s <= 0))
    stop("susie_ebnm_fn: s must be finite and positive (no exclusions).")
  s0 <- s[1]
  if (max(abs(s - s0)) > 1e-6 * s0)
    stop("susie_ebnm_fn: s must be constant across entries (use a column-wise ",
         "variance type, no missing data).")
  s2 <- s0^2

  # Effects per block: never more than the number of variants in the block.
  # (Extra effects on a tiny block would all pile onto the same variants and
  # distort the implied prior sparsity.)
  Lb <- pmin(as.integer(L), inf$nc)

  # --- collapsed data ----------------------------------------------------------
  xtx <- vector("list", B)
  for (b in multi)
    xtx[[b]] <- as.vector(crossprod(inf$X[[b]], x[inf$row_start[b]:inf$row_end[b]]))
  ns     <- length(single)
  ok_s   <- inf$sd2 > 0                                   # singleton with a usable design
  s_xtx  <- inf$sx * x[inf$srow]                          # X'x for 1x1 designs
  xh_s   <- ifelse(ok_s, s_xtx / inf$sd2, 0)
  sh2_s  <- ifelse(ok_s, s2 / inf$sd2, Inf)
  z_s    <- ifelse(ok_s, s_xtx / (sqrt(inf$sd2) * s0), 0)
  skip_s <- if (is.null(skip_z) || ns == 0) rep(FALSE, ns) else abs(z_s) < skip_z

  # --- warm start / prior ----------------------------------------------------
  has_state <- inherits(g_init, "susie_ebnm_state")
  blocks_ok <- has_state && !is.null(g_init$blocks) && length(g_init$blocks) == B &&
    all(vapply(multi, function(b)
      identical(dim(g_init$blocks[[b]]$alpha), c(Lb[b], inf$nc[b])), TRUE))

  xh0 <- c(unlist(lapply(multi, function(b) xtx[[b]] / pmax(inf$d[[b]], 1e-12))), xh_s[ok_s])
  sh0 <- c(unlist(lapply(multi, function(b) s0 / sqrt(pmax(inf$d[[b]], 1e-12)))), sqrt(sh2_s[ok_s]))
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

  blocks <- vector("list", B)                          # NULL entries = singletons
  for (b in multi) {
    if (blocks_ok) {
      blocks[[b]] <- g_init$blocks[[b]]
    } else {
      n <- inf$nc[b]; Lbb <- Lb[b]
      blocks[[b]] <- list(alpha = matrix(1 / n, Lbb, n), mu = matrix(0, Lbb, n),
                          nu = matrix(0, Lbb, n), pnz = matrix(0, Lbb, n), KL = 0,
                          counts = numeric(length(pi_k)), ell = numeric(n), V = 0,
                          skipped = FALSE)
    }
    blocks[[b]]$skipped <- FALSE
    if (!is.null(skip_z)) {
      z <- xtx[[b]] / (sqrt(pmax(inf$d[[b]], 1e-12)) * s0)
      if (max(abs(z)) < skip_z) blocks[[b]]$skipped <- TRUE
    }
  }
  sing <- list(mu = numeric(ns), nu = numeric(ns), pnz = numeric(ns), kl = numeric(ns))

  fit_pass <- function(pi_k) {
    cnt <- numeric(length(pi_k))
    for (b in multi) {
      if (blocks[[b]]$skipped) {
        n <- inf$nc[b]; Lbb <- Lb[b]
        blocks[[b]] <<- modifyList(blocks[[b]], list(
          alpha = matrix(1 / n, Lbb, n), mu = matrix(0, Lbb, n), nu = matrix(0, Lbb, n),
          pnz = matrix(0, Lbb, n), KL = -Lbb * log(pi_k[1]), ell = numeric(n), V = 0))
        blocks[[b]]$counts <<- c(Lbb, rep(0, length(pi_k) - 1))
      } else {
        blocks[[b]] <<- .block_fit(inf, b, blocks[[b]], xtx[[b]], s2, pi_k, sd_k,
                                   Lb[b], max_sweeps, alpha_tol)
      }
      cnt <- cnt + blocks[[b]]$counts
    }
    if (ns > 0) {
      sing$mu[] <<- 0; sing$nu[] <<- 0; sing$pnz[] <<- 0; sing$kl[] <<- 0
      act <- ok_s & !skip_s
      if (any(act)) {
        u <- .singleton_update(xh_s[act], sh2_s[act], pi_k, sd_k)
        sing$mu[act] <<- u$mu; sing$nu[act] <<- u$nu
        sing$pnz[act] <<- u$pnz; sing$kl[act] <<- u$kl
        cnt <- cnt + u$counts
      }
      sk <- ok_s & skip_s                              # skipped: effect fixed at 0
      if (any(sk)) {
        sing$kl[sk] <<- -log(pi_k[1])
        cnt[1] <- cnt[1] + sum(sk)
      }
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
  for (b in multi) {
    rows <- inf$row_start[b]:inf$row_end[b]
    Et[rows] <- as.vector(inf$X[[b]] %*% blocks[[b]]$ell)
    KLtot <- KLtot + blocks[[b]]$KL
    Vtot  <- Vtot + blocks[[b]]$V
  }
  if (ns > 0) {
    Et[inf$srow] <- inf$sx * sing$mu
    KLtot <- KLtot + sum(sing$kl)
    Vtot  <- Vtot + sum(inf$sd2 * pmax(sing$nu - sing$mu^2, 0))
  }
  Et2  <- Et^2 + Vtot / N
  elbo <- -0.5 * (N * log(2 * pi * s2) + (sum((x - Et)^2) + Vtot) / s2) - KLtot

  state <- structure(list(sd = sd_k, pi = pi_k, blocks = blocks, L = L, Lb = Lb,
                          single = list(idx = single, mu = sing$mu, nu = sing$nu,
                                        pnz = sing$pnz)),
                     class = "susie_ebnm_state")

  if (identical(output, "lfsr"))
    return(list(posterior = data.frame(lfsr = rep(NA_real_, N)), fitted_g = state))

  if (identical(output, "posterior_sampler")) {
    sampler <- function(nsamp) {
      out <- matrix(0, nsamp, N)
      for (b in multi) {
        bk <- blocks[[b]]; if (bk$skipped) next
        ell <- matrix(0, nsamp, inf$nc[b])
        for (m in seq_len(nrow(bk$alpha))) {
          idx <- sample.int(inf$nc[b], nsamp, replace = TRUE, prob = bk$alpha[m, ])
          sdv <- sqrt(pmax(bk$nu[m, idx] - bk$mu[m, idx]^2, 0))  # moment-matched normal
          ell[cbind(seq_len(nsamp), idx)] <- ell[cbind(seq_len(nsamp), idx)] +
            rnorm(nsamp, bk$mu[m, idx], sdv)
        }
        out[, inf$row_start[b]:inf$row_end[b]] <- ell %*% t(inf$X[[b]])
      }
      if (ns > 0) {
        sdv <- sqrt(pmax(sing$nu - sing$mu^2, 0))
        draws <- matrix(rnorm(nsamp * ns, rep(sing$mu, each = nsamp), rep(sdv, each = nsamp)),
                        nsamp, ns)
        out[, inf$srow] <- sweep(draws, 2, inf$sx, "*")
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
  out <- vector("list", length(state$blocks))
  for (b in seq_along(state$blocks)) {
    bk <- state$blocks[[b]]
    if (is.null(bk)) next
    # PIP_i = 1 - prod_m (1 - alpha_mi * P(beta_m != 0 | gamma_m = i)): unused
    # effects (mass on the point-mass-at-zero component) then add ~nothing.
    out[[b]] <- data.frame(pip = 1 - apply(1 - bk$alpha * bk$pnz, 2, prod), ell_pm = bk$ell)
  }
  sg <- state$single
  for (i in seq_along(sg$idx))                       # singleton: one effect, alpha = 1
    out[[sg$idx[i]]] <- data.frame(pip = sg$pnz[i], ell_pm = sg$mu[i])
  out
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
