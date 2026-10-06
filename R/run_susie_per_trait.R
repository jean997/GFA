run_susie_per_trait <- function(beta_hat, se, ld_list, n,
                                L = 10,
                                coverage = 0.95,
                                min_abs_corr = 0.5,
                                skip_z = 3,
                                estimate_residual_variance = FALSE,
                                ncores = 1,
                                keep_fits = FALSE,
                                verbose = TRUE) {
  beta_hat <- as.matrix(beta_hat); se <- as.matrix(se)
  stopifnot(all(dim(beta_hat) == dim(se)))
  J <- nrow(beta_hat); Tn <- ncol(beta_hat)
  trait_names <- colnames(beta_hat)
  if (is.null(trait_names)) trait_names <- paste0("T", seq_len(Tn))

  sizes  <- vapply(ld_list, nrow, 1L)
  if (sum(sizes) != J)
    stop("Sum of LD block sizes (", sum(sizes), ") != number of variants (", J, ").")
  ends   <- cumsum(sizes); starts <- ends - sizes + 1L
  nblock <- length(ld_list)

  z <- beta_hat / se
  z[!is.finite(z)] <- 0

  # Dense copies of the blocks are made lazily, once, and reused across traits.
  get_R <- local({
    cache <- vector("list", nblock)
    function(b) {
      if (is.null(cache[[b]])) cache[[b]] <<- as.matrix(ld_list[[b]])
      cache[[b]]
    }
  })

  fit_block <- function(b, t) {
    idx <- starts[b]:ends[b]
    zb  <- z[idx, t]
    if (max(abs(zb)) < skip_z || length(idx) < 2) return(NULL)
    susie_rss(z = zb, R = get_R(b), n = n, L = min(L, length(idx)),
              coverage = coverage, min_abs_corr = min_abs_corr,
              estimate_residual_variance = estimate_residual_variance)
  }

  pip  <- matrix(0, J, Tn, dimnames = list(rownames(beta_hat), trait_names))
  cs_rows <- list()
  fits <- if (keep_fits) setNames(vector("list", Tn), trait_names) else NULL

  for (t in seq_len(Tn)) {
    if (verbose) message("Trait ", trait_names[t], ": fitting ", nblock, " blocks")
    res <- if (ncores > 1) parallel::mclapply(seq_len(nblock), fit_block, t = t, mc.cores = ncores)
    else lapply(seq_len(nblock), fit_block, t = t)
    if (keep_fits) fits[[t]] <- res

    for (b in seq_len(nblock)) {
      fit <- res[[b]]
      if (is.null(fit)) next
      if (inherits(fit, "try-error")) { warning("block ", b, " failed for ", trait_names[t]); next }
      idx <- starts[b]:ends[b]
      pip[idx, t] <- fit$pip

      cs <- fit$sets$cs
      if (!is.null(cs) && length(cs) > 0) {
        pur <- fit$sets$purity
        for (k in seq_along(cs)) {
          loc  <- cs[[k]]
          glob <- idx[loc]
          lead <- glob[which.max(fit$pip[loc])]
          cs_rows[[length(cs_rows) + 1L]] <- data.frame(
            trait = trait_names[t], block = b, cs = names(cs)[k],
            size = length(loc),
            min_abs_corr = if (!is.null(pur)) pur$min.abs.corr[k] else NA_real_,
            lead = lead, lead_pip = max(fit$pip[loc]),
            stringsAsFactors = FALSE)
          cs_rows[[length(cs_rows)]]$members <- list(glob)
        }
      }
    }
    if (verbose) message("  ", sum(vapply(cs_rows, function(r) r$trait == trait_names[t], TRUE)),
                         " credible sets; max PIP = ", round(max(pip[, t]), 3))
  }

  cs_df <- if (length(cs_rows)) do.call(rbind, cs_rows) else
    data.frame(trait = character(), block = integer(), cs = character(), size = integer(),
               min_abs_corr = numeric(), lead = integer(), lead_pip = numeric())
  list(pip = pip, cs = cs_df, fits = fits)
}
